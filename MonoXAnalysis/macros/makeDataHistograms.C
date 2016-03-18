#include "makehist.h"
#include "TChain.h"

using namespace std;

void signalHiggshist(TFile* outfile,
		     int    category,
		     vector<string> observables,
		     vector<string> observables_2D,
		     double lumi                = 2.30,
		     bool   doShapeSystematics  = false,
		     string mH = "125",
		     vector<float> xs = {4.198E+04,3.925E+03,1.475E+03,9.095E+02},
		     int typeOfHiggsSignal = 0){

  if(xs.size() != 4)
    cerr<<"signalHiggshist: xs size wrong, should be 4 numbers for each mass points"<<endl;

  cout<<"Start HiggsInvisible: signalHiggshist --> "<<mH<<endl;
  
  TFile* ggHFile = NULL, *vbfHFile = NULL, *wHFile = NULL, *zHFile = NULL;
  if(typeOfHiggsSignal == 0){
    ggHFile  = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_GluGlu_HToInvisible_Mphi-"+mH+"_Mchi-0_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    vbfHFile = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_VBF_HToInvisible_Mphi-"+mH+"_Mchi-0_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    wHFile   = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_DM_ScalarWH_Mphi-"+mH+"_Mchi-10_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
    zHFile   = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_DM_ScalarZH_Mphi-"+mH+"_Mchi-10_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
  }
  else if(typeOfHiggsSignal == 1){ // fermion only
    ggHFile  = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_GluGlu_HToInvisible_Mphi-"+mH+"_Mchi-0_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    vbfHFile = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_VBF_HToInvisible_Mphi-"+mH+"_Mchi-0_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
  }
  else if(typeOfHiggsSignal == 2){  // boson only
    wHFile   = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_DM_ScalarWH_Mphi-"+mH+"_Mchi-10_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
    zHFile   = TFile::Open((baseInputTreePath+"/HiggsInvisible/sigfilter/sig_tree_DM_ScalarZH_Mphi-"+mH+"_Mchi-10_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
  }

  vector<TH1*>  ggHhist;
  vector<TH1*>  vbfHhist;
  vector<TH1*>  wHhist;
  vector<TH1*>  zHhist;
  vector<TH1*>  ggHhist_renUp, ggHhist_renDw, ggHhist_facUp, ggHhist_facDw;
  vector<TH1*>  ggHhist_bUp, ggHhist_bDw, ggHhist_metJetUp, ggHhist_metJetDw, ggHhist_metResUp, ggHhist_metResDw, ggHhist_metUncUp, ggHhist_metUncDw;
  vector<TH1*>  vbfHhist_bUp, vbfHhist_bDw, vbfHhist_metJetUp, vbfHhist_metJetDw, vbfHhist_metResUp, vbfHhist_metResDw, vbfHhist_metUncUp, vbfHhist_metUncDw;
  vector<TH1*>  wHhist_bUp, wHhist_bDw, wHhist_metJetUp, wHhist_metJetDw, wHhist_metResUp, wHhist_metResDw, wHhist_metUncUp, wHhist_metUncDw;
  vector<TH1*>  zHhist_bUp, zHhist_bDw, zHhist_metJetUp, zHhist_metJetDw, zHhist_metResUp, zHhist_metResDw, zHhist_metUncUp, zHhist_metUncDw;
  vector<float> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(typeOfHiggsSignal <= 1){
      TH1F* ggHhist_temp  = new TH1F(("ggHhist_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
      TH1F* vbfHhist_temp = new TH1F(("vbfHhist_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      ggHhist.push_back(dynamic_cast<TH1*>(ggHhist_temp));
      vbfHhist.push_back(dynamic_cast<TH1*>(vbfHhist_temp));
    }
    if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
      TH1F* wHhist_temp   = new TH1F(("wHhist_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
      TH1F* zHhist_temp   = new TH1F(("zHhist_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
      wHhist.push_back(dynamic_cast<TH1*>(wHhist_temp));
      zHhist.push_back(dynamic_cast<TH1*>(zHhist_temp));
    }

    if(doShapeSystematics){

      if(typeOfHiggsSignal <= 1){
	TH1F* ggHhist_renUp_temp   = new TH1F(("ggHhist_renUp_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* ggHhist_renDw_temp   = new TH1F(("ggHhist_renDw_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* ggHhist_facUp_temp   = new TH1F(("ggHhist_facUp_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* ggHhist_facDw_temp   = new TH1F(("ggHhist_facDw_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);      
	ggHhist_renUp.push_back(dynamic_cast<TH1*>(ggHhist_renUp_temp));
	ggHhist_renDw.push_back(dynamic_cast<TH1*>(ggHhist_renDw_temp));
	ggHhist_facUp.push_back(dynamic_cast<TH1*>(ggHhist_facUp_temp));
	ggHhist_facDw.push_back(dynamic_cast<TH1*>(ggHhist_facDw_temp));

	TH1F* ggHhist_bUp_temp  = new TH1F(("ggHhist_bUp_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_bUp_temp = new TH1F(("vbfHhist_bUp_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_bUp.push_back(dynamic_cast<TH1*>(ggHhist_bUp_temp));
	vbfHhist_bUp.push_back(dynamic_cast<TH1*>(vbfHhist_bUp_temp));
      
	TH1F* ggHhist_bDw_temp  = new TH1F(("ggHhist_bDw_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_bDw_temp = new TH1F(("vbfHhist_bDw_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_bDw.push_back(dynamic_cast<TH1*>(ggHhist_bDw_temp));
	vbfHhist_bDw.push_back(dynamic_cast<TH1*>(vbfHhist_bDw_temp));

	TH1F* ggHhist_metJetUp_temp  = new TH1F(("ggHhist_metJetUp_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_metJetUp_temp = new TH1F(("vbfHhist_metJetUp_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_metJetUp.push_back(dynamic_cast<TH1*>(ggHhist_metJetUp_temp));
	vbfHhist_metJetUp.push_back(dynamic_cast<TH1*>(vbfHhist_metJetUp_temp));

	TH1F* ggHhist_metJetDw_temp  = new TH1F(("ggHhist_metJetDw_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_metJetDw_temp = new TH1F(("vbfHhist_metJetDw_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_metJetDw.push_back(dynamic_cast<TH1*>(ggHhist_metJetDw_temp));
	vbfHhist_metJetDw.push_back(dynamic_cast<TH1*>(vbfHhist_metJetDw_temp));

	TH1F* ggHhist_metResUp_temp  = new TH1F(("ggHhist_metResUp_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_metResUp_temp = new TH1F(("vbfHhist_metResUp_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_metResUp.push_back(dynamic_cast<TH1*>(ggHhist_metResUp_temp));
	vbfHhist_metResUp.push_back(dynamic_cast<TH1*>(vbfHhist_metResUp_temp));

	TH1F* ggHhist_metResDw_temp  = new TH1F(("ggHhist_metResDw_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_metResDw_temp = new TH1F(("vbfHhist_metResDw_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_metResDw.push_back(dynamic_cast<TH1*>(ggHhist_metResDw_temp));
	vbfHhist_metResDw.push_back(dynamic_cast<TH1*>(vbfHhist_metResDw_temp));

	TH1F* ggHhist_metUncUp_temp  = new TH1F(("ggHhist_metUncUp_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_metUncUp_temp = new TH1F(("vbfHhist_metUncUp_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_metUncUp.push_back(dynamic_cast<TH1*>(ggHhist_metUncUp_temp));
	vbfHhist_metUncUp.push_back(dynamic_cast<TH1*>(vbfHhist_metUncUp_temp));

	TH1F* ggHhist_metUncDw_temp  = new TH1F(("ggHhist_metUncDw_"+mH+"_"+obs).c_str(), "",  int(bins.size()-1), &bins[0]);
	TH1F* vbfHhist_metUncDw_temp = new TH1F(("vbfHhist_metUncDw_"+mH+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	ggHhist_metUncDw.push_back(dynamic_cast<TH1*>(ggHhist_metUncDw_temp));
	vbfHhist_metUncDw.push_back(dynamic_cast<TH1*>(vbfHhist_metUncDw_temp));
      }
      else if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
	
	TH1F* wHhist_bUp_temp   = new TH1F(("wHhist_bUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_bUp_temp   = new TH1F(("zHhist_bUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_bUp.push_back(dynamic_cast<TH1*>(wHhist_bUp_temp));
	zHhist_bUp.push_back(dynamic_cast<TH1*>(zHhist_bUp_temp));

	TH1F* wHhist_bDw_temp   = new TH1F(("wHhist_bDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_bDw_temp   = new TH1F(("zHhist_bDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_bDw.push_back(dynamic_cast<TH1*>(wHhist_bDw_temp));
	zHhist_bDw.push_back(dynamic_cast<TH1*>(zHhist_bDw_temp));
	
	TH1F* wHhist_metJetUp_temp   = new TH1F(("wHhist_metJetUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_metJetUp_temp   = new TH1F(("zHhist_metJetUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	zHhist_metJetUp.push_back(dynamic_cast<TH1*>(zHhist_metJetUp_temp));
	wHhist_metJetUp.push_back(dynamic_cast<TH1*>(wHhist_metJetUp_temp));
	
	TH1F* wHhist_metJetDw_temp   = new TH1F(("wHhist_metJetDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_metJetDw_temp   = new TH1F(("zHhist_metJetDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_metJetDw.push_back(dynamic_cast<TH1*>(wHhist_metJetDw_temp));
	zHhist_metJetDw.push_back(dynamic_cast<TH1*>(zHhist_metJetDw_temp));
	
	TH1F* wHhist_metResUp_temp   = new TH1F(("wHhist_metResUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_metResUp_temp   = new TH1F(("zHhist_metResUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_metResUp.push_back(dynamic_cast<TH1*>(wHhist_metResUp_temp));
	zHhist_metResUp.push_back(dynamic_cast<TH1*>(zHhist_metResUp_temp));
	
	TH1F* wHhist_metResDw_temp   = new TH1F(("wHhist_metResDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_metResDw_temp   = new TH1F(("zHhist_metResDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_metResDw.push_back(dynamic_cast<TH1*>(wHhist_metResDw_temp));
	zHhist_metResDw.push_back(dynamic_cast<TH1*>(zHhist_metResDw_temp));
	
	TH1F* wHhist_metUncUp_temp   = new TH1F(("wHhist_metUncUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_metUncUp_temp   = new TH1F(("zHhist_metUncUp_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_metUncUp.push_back(dynamic_cast<TH1*>(wHhist_metUncUp_temp));
	zHhist_metUncUp.push_back(dynamic_cast<TH1*>(zHhist_metUncUp_temp));
	
	TH1F* wHhist_metUncDw_temp   = new TH1F(("wHhist_metUncDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	TH1F* zHhist_metUncDw_temp   = new TH1F(("zHhist_metUncDw_"+mH+"_"+obs).c_str(), "",   int(bins.size()-1), &bins[0]);
	wHhist_metUncDw.push_back(dynamic_cast<TH1*>(wHhist_metUncDw_temp));
	zHhist_metUncDw.push_back(dynamic_cast<TH1*>(zHhist_metUncDw_temp));
	
      }  
    }
  }

  vector<TH2*> ggHhist_2D; 
  vector<TH2*> vbfHhist_2D;
  vector<TH2*> wHhist_2D;
  vector<TH2*> zHhist_2D;

  vector<TH2*> ggHhist_renUp_2D, ggHhist_renDw_2D, ggHhist_facUp_2D, ggHhist_facDw_2D;
  vector<TH2*> ggHhist_bUp_2D, ggHhist_bDw_2D, ggHhist_metJetUp_2D, ggHhist_metJetDw_2D, ggHhist_metResUp_2D, ggHhist_metResDw_2D, ggHhist_metUncUp_2D, ggHhist_metUncDw_2D;
  vector<TH2*> vbfHhist_bUp_2D, vbfHhist_bDw_2D, vbfHhist_metJetUp_2D, vbfHhist_metJetDw_2D, vbfHhist_metResUp_2D, vbfHhist_metResDw_2D, vbfHhist_metUncUp_2D, vbfHhist_metUncDw_2D;
  vector<TH2*> wHhist_bUp_2D, wHhist_bDw_2D, wHhist_metJetUp_2D, wHhist_metJetDw_2D, wHhist_metResUp_2D, wHhist_metResDw_2D, wHhist_metUncUp_2D, wHhist_metUncDw_2D;
  vector<TH2*> zHhist_bUp_2D, zHhist_bDw_2D, zHhist_metJetUp_2D, zHhist_metJetDw_2D, zHhist_metResUp_2D, zHhist_metResDw_2D, zHhist_metUncUp_2D, zHhist_metUncDw_2D;

  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    if(typeOfHiggsSignal <= 1){
      TH2F* ggHhist_temp  = new TH2F(("ggHhist_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vbfHhist_temp = new TH2F(("vbfHhist_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      ggHhist_2D.push_back(dynamic_cast<TH2*>(ggHhist_temp));
      vbfHhist_2D.push_back(dynamic_cast<TH2*>(vbfHhist_temp));
    }
    else if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
      TH2F* wHhist_temp   = new TH2F(("wHhist_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* zHhist_temp   = new TH2F(("zHhist_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      wHhist_2D.push_back(dynamic_cast<TH2*>(wHhist_temp));
      zHhist_2D.push_back(dynamic_cast<TH2*>(zHhist_temp));
    }
    
    if(doShapeSystematics){

      if(typeOfHiggsSignal <= 1){

	TH2F* ggHhist_renUp_temp   = new TH2F(("ggHhist_renUp_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* ggHhist_renDw_temp   = new TH2F(("ggHhist_renDw_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* ggHhist_facUp_temp   = new TH2F(("ggHhist_facUp_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* ggHhist_facDw_temp   = new TH2F(("ggHhist_facDw_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_renUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_renUp_temp));
	ggHhist_renDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_renDw_temp));
	ggHhist_facUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_facUp_temp));
	ggHhist_facDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_facDw_temp));
      
	TH2F* ggHhist_bUp_temp  = new TH2F(("ggHhist_bUp_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_bUp_temp = new TH2F(("vbfHhist_bUp_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_bUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_bUp_temp));
	vbfHhist_bUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_bUp_temp));

	TH2F* ggHhist_bDw_temp  = new TH2F(("ggHhist_bDw_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_bDw_temp = new TH2F(("vbfHhist_bDw_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_bDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_bDw_temp));
	vbfHhist_bDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_bDw_temp));

	TH2F* ggHhist_metJetUp_temp  = new TH2F(("ggHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_metJetUp_temp = new TH2F(("vbfHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_metJetUp_temp));
	vbfHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metJetUp_temp));

	TH2F* ggHhist_metJetDw_temp  = new TH2F(("ggHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_metJetDw_temp = new TH2F(("vbfHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_metJetDw_temp));
	vbfHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metJetDw_temp));

	TH2F* ggHhist_metResUp_temp  = new TH2F(("ggHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_metResUp_temp = new TH2F(("vbfHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_metResUp_temp));
	vbfHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metResUp_temp));
	
	TH2F* ggHhist_metResDw_temp  = new TH2F(("ggHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_metResDw_temp = new TH2F(("vbfHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_metResDw_temp));
	vbfHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metResDw_temp));
	
	TH2F* ggHhist_metUncUp_temp  = new TH2F(("ggHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_metUncUp_temp = new TH2F(("vbfHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(ggHhist_metUncUp_temp));
	vbfHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metUncUp_temp));

	TH2F* ggHhist_metUncDw_temp  = new TH2F(("ggHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(), "",  int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* vbfHhist_metUncDw_temp = new TH2F(("vbfHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	ggHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(ggHhist_metUncDw_temp));
	vbfHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vbfHhist_metUncDw_temp));
      }
      else if(typeOfHiggsSignal == 0 or typeOfHiggsSignal>= 2){

	TH2F* wHhist_bUp_temp   = new TH2F(("wHhist_bUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_bUp_temp   = new TH2F(("zHhist_bUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	wHhist_bUp_2D.push_back(dynamic_cast<TH2*>(wHhist_bUp_temp));
	zHhist_bUp_2D.push_back(dynamic_cast<TH2*>(zHhist_bUp_temp));
	
	TH2F* wHhist_bDw_temp   = new TH2F(("wHhist_bDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_bDw_temp   = new TH2F(("zHhist_bDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	wHhist_bDw_2D.push_back(dynamic_cast<TH2*>(wHhist_bDw_temp));
	zHhist_bDw_2D.push_back(dynamic_cast<TH2*>(zHhist_bDw_temp));

	TH2F* wHhist_metJetUp_temp   = new TH2F(("wHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_metJetUp_temp   = new TH2F(("zHhist_metJetUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	zHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(zHhist_metJetUp_temp));
	wHhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(wHhist_metJetUp_temp));
	
	TH2F* wHhist_metJetDw_temp   = new TH2F(("wHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_metJetDw_temp   = new TH2F(("zHhist_metJetDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	wHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(wHhist_metJetDw_temp));
	zHhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(zHhist_metJetDw_temp));
	
	TH2F* wHhist_metResUp_temp   = new TH2F(("wHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_metResUp_temp   = new TH2F(("zHhist_metResUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	wHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(wHhist_metResUp_temp));
	zHhist_metResUp_2D.push_back(dynamic_cast<TH2*>(zHhist_metResUp_temp));

	TH2F* wHhist_metResDw_temp   = new TH2F(("wHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_metResDw_temp   = new TH2F(("zHhist_metResDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	wHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(wHhist_metResDw_temp));
	zHhist_metResDw_2D.push_back(dynamic_cast<TH2*>(zHhist_metResDw_temp));
	
	TH2F* wHhist_metUncUp_temp   = new TH2F(("wHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_metUncUp_temp   = new TH2F(("zHhist_metUncUp_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	wHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(wHhist_metUncUp_temp));
	zHhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(zHhist_metUncUp_temp));
	
	TH2F* wHhist_metUncDw_temp   = new TH2F(("wHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* zHhist_metUncDw_temp   = new TH2F(("zHhist_metUncDw_"+mH+"_"+obs+"_2D").c_str(), "",   int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

	wHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(wHhist_metUncDw_temp));
	zHhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(zHhist_metUncDw_temp));
      }
    }  
  }
  
  // start running on signal samples                                                                                                                                 
  vector<TH1*> ehists;
  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  TTree* ggHTree = NULL, *vbfHTree = NULL, *wHTree = NULL, *zHTree = NULL;
  if(ggHFile)  ggHTree =  (TTree*) ggHFile->Get("tree/tree");
  if(vbfHFile) vbfHTree = (TTree*) vbfHFile->Get("tree/tree");
  if(wHFile)   wHTree =  (TTree*) wHFile->Get("tree/tree");
  if(zHFile)   zHTree =  (TTree*) zHFile->Get("tree/tree");

  //Take histograms for Higgs
  TFile* higgsPT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/HiggsPT/hpt_ggH_125_13TeV_POWHEG_weights.root");

  TH1* renUp = (TH1*) higgsPT->Get("weight_pT_rup");
  TH1* renDw = (TH1*) higgsPT->Get("weight_pT_rdn");
  TH1* facUp = (TH1*) higgsPT->Get("weight_pT_fup");
  TH1* facDw = (TH1*) higgsPT->Get("weight_pT_fup");
  for(int iBinX = 1; iBinX <= facDw->GetNbinsX(); iBinX++)
    facDw->SetBinContent(iBinX,1-facDw->GetBinContent(iBinX));

  // higgs re-weight
  if(typeOfHiggsSignal <=1){
    makehist4(ggHTree,ggHhist,ggHhist_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(0));
    makehist4(vbfHTree,vbfHhist,vbfHhist_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(1));
  }
  else if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >=2){
    makehist4(wHTree,wHhist,wHhist_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(2));
    makehist4(zHTree,zHhist,zHhist_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(3));
  }

  if(doShapeSystematics){
 
    if(typeOfHiggsSignal <= 1){
      makehist4(ggHTree,ggHhist_renUp,ggHhist_renUp_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(0),renUp);
      makehist4(ggHTree,ggHhist_renDw,ggHhist_renDw_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(0),renDw);
      makehist4(ggHTree,ggHhist_facUp,ggHhist_facUp_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(0),facUp);
      makehist4(ggHTree,ggHhist_facDw,ggHhist_facDw_2D,true,0,category,false,1.00,lumi,0,ehists,"",false,true,0,true,xs.at(0),facDw);
    
      makehist4(ggHTree,ggHhist_bUp,ggHhist_bUp_2D,true,0,category,false,1.00,lumi,0,ehists,"btagUp",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_bUp,vbfHhist_bUp_2D,true,0,category,false,1.00,lumi,0,ehists,"btagUp",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_bDw,ggHhist_bDw_2D,true,0,category,false,1.00,lumi,0,ehists,"btagDown",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_bDw,vbfHhist_bDw_2D,true,0,category,false,1.00,lumi,0,ehists,"btagDown",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_metJetUp,ggHhist_metJetUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jesUp",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metJetUp,vbfHhist_metJetUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jesUp",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_metJetDw,ggHhist_metJetDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jesDw",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metJetDw,vbfHhist_metJetDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jesDw",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_metResUp,ggHhist_metResUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jerUp",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metResUp,vbfHhist_metResUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jerUp",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_metResDw,ggHhist_metResDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jerDw",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metResDw,vbfHhist_metResDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jerDw",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_metUncUp,ggHhist_metUncUp_2D,true,0,category,false,1.00,lumi,0,ehists,"uncUp",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metUncUp,vbfHhist_metUncUp_2D,true,0,category,false,1.00,lumi,0,ehists,"uncUp",false,true,0,true,xs.at(1));
      makehist4(ggHTree,ggHhist_metUncDw,ggHhist_metUncDw_2D,true,0,category,false,1.00,lumi,0,ehists,"uncDw",false,true,0,true,xs.at(0));
      makehist4(vbfHTree,vbfHhist_metUncDw,vbfHhist_metUncDw_2D,true,0,category,false,1.00,lumi,0,ehists,"uncDw",false,true,0,true,xs.at(1));
    }
    else if(typeOfHiggsSignal == 0 or typeOfHiggsSignal >= 2){

      makehist4(wHTree,wHhist_bUp,wHhist_bUp_2D,true,0,category,false,1.00,lumi,0,ehists,"btagUp",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_bUp,zHhist_bUp_2D,true,0,category,false,1.00,lumi,0,ehists,"btagUp",false,true,0,true,xs.at(3));      
      makehist4(wHTree,wHhist_bDw,wHhist_bDw_2D,true,0,category,false,1.00,lumi,0,ehists,"btagDown",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_bDw,zHhist_bDw_2D,true,0,category,false,1.00,lumi,0,ehists,"btagDown",false,true,0,true,xs.at(3));   
      makehist4(wHTree,wHhist_metJetUp,wHhist_metJetUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jesUp",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_metJetUp,zHhist_metJetUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jesUp",false,true,0,true,xs.at(3));    
      makehist4(wHTree,wHhist_metJetDw,wHhist_metJetDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jesDw",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_metJetDw,zHhist_metJetDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jesDw",false,true,0,true,xs.at(3));    
      makehist4(wHTree,wHhist_metResUp,wHhist_metResUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jerUp",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_metResUp,zHhist_metResUp_2D,true,0,category,false,1.00,lumi,0,ehists,"jerUp",false,true,0,true,xs.at(3));
      makehist4(wHTree,wHhist_metResDw,wHhist_metResDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jerDw",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_metResDw,zHhist_metResDw_2D,true,0,category,false,1.00,lumi,0,ehists,"jerDw",false,true,0,true,xs.at(3));
      makehist4(wHTree,wHhist_metUncUp,wHhist_metUncUp_2D,true,0,category,false,1.00,lumi,0,ehists,"uncUp",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_metUncUp,zHhist_metUncUp_2D,true,0,category,false,1.00,lumi,0,ehists,"uncUp",false,true,0,true,xs.at(3));
      makehist4(wHTree,wHhist_metUncDw,wHhist_metUncDw_2D,true,0,category,false,1.00,lumi,0,ehists,"uncDw",false,true,0,true,xs.at(2));
      makehist4(zHTree,zHhist_metUncDw,zHhist_metUncDw_2D,true,0,category,false,1.00,lumi,0,ehists,"uncDw",false,true,0,true,xs.at(3));
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

  }

  // store tempaltes in the output file                                                                                                                                     
  outfile->cd();
  for(auto hist : ggHhist)  hist->Write();
  for(auto hist : vbfHhist) hist->Write();
  for(auto hist : wHhist)   hist->Write();
  for(auto hist : zHhist)   hist->Write();

  for(auto hist : ggHhist_2D)  hist->Write();
  for(auto hist : vbfHhist_2D) hist->Write();
  for(auto hist : wHhist_2D)   hist->Write();
  for(auto hist : zHhist_2D)   hist->Write();

  for(auto hist : ggHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
  for(auto hist : vbfHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
  for(auto hist : wHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
  for(auto hist : zHhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

  if(doShapeSystematics){

    for(auto hist : ggHhist_renUp) hist->Write();
    for(auto hist : ggHhist_renDw) hist->Write();
    for(auto hist : ggHhist_facUp) hist->Write();
    for(auto hist : ggHhist_facDw) hist->Write();

    for(auto hist : ggHhist_renUp_2D) hist->Write();
    for(auto hist : ggHhist_renDw_2D) hist->Write();
    for(auto hist : ggHhist_facUp_2D) hist->Write();
    for(auto hist : ggHhist_facDw_2D) hist->Write();

    for(auto hist : ggHhist_renUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_renDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_facUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : ggHhist_facDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_bUp) hist->Write();
    for(auto hist : vbfHhist_bUp) hist->Write();
    for(auto hist : wHhist_bUp) hist->Write();
    for(auto hist : zHhist_bUp) hist->Write();

    for(auto hist : ggHhist_bUp_2D) hist->Write();
    for(auto hist : vbfHhist_bUp_2D) hist->Write();
    for(auto hist : wHhist_bUp_2D) hist->Write();
    for(auto hist : zHhist_bUp_2D) hist->Write();

    for(auto hist : ggHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    //
    for(auto hist : ggHhist_bDw) hist->Write();
    for(auto hist : vbfHhist_bDw) hist->Write();
    for(auto hist : wHhist_bDw) hist->Write();
    for(auto hist : zHhist_bDw) hist->Write();

    for(auto hist : ggHhist_bDw_2D) hist->Write();
    for(auto hist : vbfHhist_bDw_2D) hist->Write();
    for(auto hist : wHhist_bDw_2D) hist->Write();
    for(auto hist : zHhist_bDw_2D) hist->Write();

    for(auto hist : ggHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_metJetUp) hist->Write();
    for(auto hist : vbfHhist_metJetUp) hist->Write();
    for(auto hist : wHhist_metJetUp) hist->Write();
    for(auto hist : zHhist_metJetUp) hist->Write();

    for(auto hist : ggHhist_metJetUp_2D) hist->Write();
    for(auto hist : vbfHhist_metJetUp_2D) hist->Write();
    for(auto hist : wHhist_metJetUp_2D) hist->Write();
    for(auto hist : zHhist_metJetUp_2D) hist->Write();

    for(auto hist : ggHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist :vbfHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist : wHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_metJetDw) hist->Write();
    for(auto hist : vbfHhist_metJetDw) hist->Write();
    for(auto hist : wHhist_metJetDw) hist->Write();
    for(auto hist : zHhist_metJetDw) hist->Write();

    for(auto hist : ggHhist_metJetDw_2D) hist->Write();
    for(auto hist : vbfHhist_metJetDw_2D) hist->Write();
    for(auto hist : wHhist_metJetDw_2D) hist->Write();
    for(auto hist : zHhist_metJetDw_2D) hist->Write();

    for(auto hist : ggHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_metResUp) hist->Write();
    for(auto hist : vbfHhist_metResUp) hist->Write();
    for(auto hist : wHhist_metResUp) hist->Write();
    for(auto hist : zHhist_metResUp) hist->Write();

    for(auto hist : ggHhist_metResUp_2D) hist->Write();
    for(auto hist : vbfHhist_metResUp_2D) hist->Write();
    for(auto hist : wHhist_metResUp_2D) hist->Write();
    for(auto hist : zHhist_metResUp_2D) hist->Write();

    for(auto hist : ggHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_metResDw) hist->Write();
    for(auto hist : vbfHhist_metResDw) hist->Write();
    for(auto hist : wHhist_metResDw) hist->Write();
    for(auto hist : zHhist_metResDw) hist->Write();

    for(auto hist : ggHhist_metResDw_2D) hist->Write();
    for(auto hist : vbfHhist_metResDw_2D) hist->Write();
    for(auto hist : wHhist_metResDw_2D) hist->Write();
    for(auto hist : zHhist_metResDw_2D) hist->Write();

    for(auto hist : ggHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_metUncUp) hist->Write();
    for(auto hist : vbfHhist_metUncUp) hist->Write();
    for(auto hist : wHhist_metUncUp) hist->Write();
    for(auto hist : zHhist_metUncUp) hist->Write();

    for(auto hist : ggHhist_metUncUp_2D) hist->Write();
    for(auto hist : vbfHhist_metUncUp_2D) hist->Write();
    for(auto hist : wHhist_metUncUp_2D) hist->Write();
    for(auto hist : zHhist_metUncUp_2D) hist->Write();

    for(auto hist : ggHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :vbfHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : wHhist_metUncUp_2D){  TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist : zHhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist : ggHhist_metUncDw) hist->Write();
    for(auto hist : vbfHhist_metUncDw) hist->Write();
    for(auto hist : wHhist_metUncDw) hist->Write();
    for(auto hist : zHhist_metUncDw) hist->Write();

    for(auto hist : ggHhist_metUncDw_2D) hist->Write();
    for(auto hist : vbfHhist_metUncDw_2D) hist->Write();
    for(auto hist : wHhist_metUncDw_2D) hist->Write();
    for(auto hist : zHhist_metUncDw_2D) hist->Write();

    for(auto hist : ggHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist :vbfHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist : wHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}
    for(auto hist : zHhist_metUncDw_2D){TH1* temp = unroll2DHistograms(hist);temp->Write();}
  }

  ggHhist.clear();
  vbfHhist.clear();
  wHhist.clear();
  zHhist.clear();
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

  if(ggHFile) ggHFile->Close();
  if(vbfHFile) vbfHFile->Close();
  if(wHFile) wHFile->Close();
  if(zHFile) zHFile->Close();
  if(higgsPT) higgsPT->Close();

  cout << "Templates for the Higgs invisible ..." << endl;
}


// Build templates for the signal region                                                                                                                                       
void signalmchist(TFile* outfile,
		  int category,
		  vector<string> observables,
		  vector<string> observables_2D,
		  vector<signalSample> massPoint,
		  string interaction       = "Vector",
		  double lumi              = 2.24,
		  bool doShapeSystematics  = false){

  
  vector<TFile*> monoJfile;
  vector<TFile*> monoWfile;
  vector<TFile*> monoZfile;

  for(auto iPoint : massPoint){
    if(iPoint.interaction == interaction and interaction == "Vector"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMV_Vector/sigfilter/sig_tree_DMV_NNPDF30_Vector_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
      monoZfile .push_back(TFile::Open((baseInputTreePath+"MonoZ_Vector/sigfilter/sig_tree_VectorMonoZ_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
    }
    else if(iPoint.interaction == interaction and interaction == "Axial"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMV_Axial/sigfilter/sig_tree_DMV_NNPDF30_Axial_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
      monoZfile .push_back(TFile::Open((baseInputTreePath+"MonoZ_Axial/sigfilter/sig_tree_AxialMonoZ_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str()));
    }
    else if(iPoint.interaction == interaction and interaction == "Scalar"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMS_Scalar/sigfilter/sig_tree_DMS_NNPDF30_Scalar_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-"+iPoint.mediatorMass+
					"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str()));      
    }
    else if(iPoint.interaction == interaction and interaction == "Pseudoscalar"){
      monoJfile .push_back(TFile::Open((baseInputTreePath+"DMS_Pseudoscalar/sigfilter/sig_tree_DMS_NNPDF30_Pseudoscalar_Mphi-"+
					iPoint.mediatorMass+"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str()));
      monoWfile .push_back(TFile::Open((baseInputTreePath+"MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-"+
					iPoint.mediatorMass+"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str()));
      monoZfile .push_back(TFile::Open((baseInputTreePath+"MonoZ_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarZH_Mphi-"+
					iPoint.mediatorMass+"_Mchi-"+iPoint.dmMass+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str()));
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
  vector<float> bins;
  for(auto iPoint : massPoint){
    if(iPoint.interaction != interaction) continue;
    
    for(auto obs : observables){
      
      bins = selectBinning(obs,category);
      if(bins.empty())
        cout<<"No binning for this observable --> please define it"<<endl;
 
      TH1F* monoJhist_temp = new TH1F(("monoJhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoWhist_temp = new TH1F(("monoWhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoZhist_temp = new TH1F(("monoZhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      monoJhist.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_temp));
      monoWhist.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_temp));
      monoZhist.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_temp));

      if(doShapeSystematics){

	TH1F* monoJhist_bUp_temp = new TH1F(("monoJhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_bUp_temp = new TH1F(("monoWhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_bUp_temp = new TH1F(("monoZhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_bUp_temp));
	monoWhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_bUp_temp));
	monoZhist_bUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_bUp_temp));

	TH1F* monoJhist_bDw_temp = new TH1F(("monoJhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_bDw_temp = new TH1F(("monoWhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_bDw_temp = new TH1F(("monoZhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_bDw_temp));
	monoWhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_bDw_temp));
	monoZhist_bDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_bDw_temp));


	TH1F* monoJhist_metJetUp_temp = new TH1F(("monoJhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metJetUp_temp = new TH1F(("monoWhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metJetUp_temp = new TH1F(("monoZhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metJetUp_temp));
	monoWhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metJetUp_temp));
	monoZhist_metJetUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metJetUp_temp));

	TH1F* monoJhist_metJetDw_temp = new TH1F(("monoJhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metJetDw_temp = new TH1F(("monoWhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metJetDw_temp = new TH1F(("monoZhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metJetDw_temp));
	monoWhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metJetDw_temp));
	monoZhist_metJetDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metJetDw_temp));

	TH1F* monoJhist_metResUp_temp = new TH1F(("monoJhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metResUp_temp = new TH1F(("monoWhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metResUp_temp = new TH1F(("monoZhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metResUp_temp));
	monoWhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metResUp_temp));
	monoZhist_metResUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metResUp_temp));

	TH1F* monoJhist_metResDw_temp = new TH1F(("monoJhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metResDw_temp = new TH1F(("monoWhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metResDw_temp = new TH1F(("monoZhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metResDw_temp));
	monoWhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metResDw_temp));
	monoZhist_metResDw.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metResDw_temp));

	TH1F* monoJhist_metUncUp_temp = new TH1F(("monoJhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metUncUp_temp = new TH1F(("monoWhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metUncUp_temp = new TH1F(("monoZhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	monoJhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoJhist_metUncUp_temp));
	monoWhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoWhist_metUncUp_temp));
	monoZhist_metUncUp.at(imass).push_back(dynamic_cast<TH1*>(monoZhist_metUncUp_temp));

	TH1F* monoJhist_metUncDw_temp = new TH1F(("monoJhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoWhist_metUncDw_temp = new TH1F(("monoWhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
	TH1F* monoZhist_metUncDw_temp = new TH1F(("monoZhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
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

      TH2F* monoJhist_temp = new TH2F(("monoJhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* monoWhist_temp = new TH2F(("monoWhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* monoZhist_temp = new TH2F(("monoZhist_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      monoJhist_2D.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_temp));
      monoWhist_2D.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_temp));
      monoZhist_2D.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_temp));

      if(doShapeSystematics){

	TH2F* monoJhist_bUp_temp = new TH2F(("monoJhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_bUp_temp = new TH2F(("monoWhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_bUp_temp = new TH2F(("monoZhist_bUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_bUp.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_bUp_temp));
	monoWhist_bUp.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_bUp_temp));
	monoZhist_bUp.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_bUp_temp));

	TH2F* monoJhist_bDw_temp = new TH2F(("monoJhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_bDw_temp = new TH2F(("monoWhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_bDw_temp = new TH2F(("monoZhist_bDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_bDw.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_bDw_temp));
	monoWhist_bDw.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_bDw_temp));
	monoZhist_bDw.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_bDw_temp));


	TH2F* monoJhist_metJetUp_temp = new TH2F(("monoJhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_metJetUp_temp = new TH2F(("monoWhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_metJetUp_temp = new TH2F(("monoZhist_metJetUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_metJetUp.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metJetUp_temp));
	monoWhist_metJetUp.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metJetUp_temp));
	monoZhist_metJetUp.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metJetUp_temp));

	TH2F* monoJhist_metJetDw_temp = new TH2F(("monoJhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_metJetDw_temp = new TH2F(("monoWhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_metJetDw_temp = new TH2F(("monoZhist_metJetDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_metJetDw.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metJetDw_temp));
	monoWhist_metJetDw.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metJetDw_temp));
	monoZhist_metJetDw.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metJetDw_temp));

	TH2F* monoJhist_metResUp_temp = new TH2F(("monoJhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_metResUp_temp = new TH2F(("monoWhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_metResUp_temp = new TH2F(("monoZhist_metResUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_metResUp.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metResUp_temp));
	monoWhist_metResUp.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metResUp_temp));
	monoZhist_metResUp.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metResUp_temp));

	TH2F* monoJhist_metResDw_temp = new TH2F(("monoJhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_metResDw_temp = new TH2F(("monoWhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_metResDw_temp = new TH2F(("monoZhist_metResDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_metResDw.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metResDw_temp));
	monoWhist_metResDw.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metResDw_temp));
	monoZhist_metResDw.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metResDw_temp));

	TH2F* monoJhist_metUncUp_temp = new TH2F(("monoJhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_metUncUp_temp = new TH2F(("monoWhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_metUncUp_temp = new TH2F(("monoZhist_metUncUp_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_metUncUp.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metUncUp_temp));
	monoWhist_metUncUp.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metUncUp_temp));
	monoZhist_metUncUp.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metUncUp_temp));

	TH2F* monoJhist_metUncDw_temp = new TH2F(("monoJhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoWhist_metUncDw_temp = new TH2F(("monoWhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	TH2F* monoZhist_metUncDw_temp = new TH2F(("monoZhist_metUncDw_"+interaction+"_"+iPoint.mediatorMass+"_"+iPoint.dmMass+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
	monoJhist_metUncDw.at(imass).push_back(dynamic_cast<TH2*>(monoJhist_metUncDw_temp));
	monoWhist_metUncDw.at(imass).push_back(dynamic_cast<TH2*>(monoWhist_metUncDw_temp));
	monoZhist_metUncDw.at(imass).push_back(dynamic_cast<TH2*>(monoZhist_metUncDw_temp));

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
  if(category == 2 or category == 3)
    isWJet = true;

  int itree = 0;
  for(auto tree : monoJtree){
    // signals                                                                                                                                                                 
    if(tree){
      cout<<"signal region analysis --> Signal monoJet "<<endl;
      makehist4(tree, monoJhist.at(itree),  monoJhist_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "",  false, true);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoJ sys "<<endl;
	makehist4(tree, monoJhist_bUp.at(itree),  monoJhist_bUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "btagUp",  false, true);
	makehist4(tree, monoJhist_bDw.at(itree),  monoJhist_bDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "btagDown",  false, true);
	makehist4(tree, monoJhist_metJetUp.at(itree),  monoJhist_metJetUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jesUp", false, true);
	makehist4(tree, monoJhist_metJetDw.at(itree),  monoJhist_metJetDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jesDw", false, true);
	makehist4(tree, monoJhist_metResUp.at(itree),  monoJhist_metResUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerUp", false, true);
	makehist4(tree, monoJhist_metResDw.at(itree),  monoJhist_metResDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "jerDw", false, true);
	makehist4(tree, monoJhist_metUncUp.at(itree),  monoJhist_metUncUp_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncUp", false, true);
	makehist4(tree, monoJhist_metUncDw.at(itree),  monoJhist_metUncDw_2D.at(itree),  true, 0, category, false, 1.00, lumi, 0, ehists, "uncDw", false, true);
      }
    }
    itree++;
  }

  itree = 0;
  for(auto tree : monoWtree){
    if(tree){
      cout<<"signal region analysis --> Signal monoW "<<endl;
      makehist4(tree, monoWhist.at(itree),  monoWhist_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "", false, true);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoW sys "<<endl;
	makehist4(tree, monoWhist_bUp.at(itree),  monoWhist_bUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "btagUp",  false, true);
	makehist4(tree, monoWhist_bDw.at(itree),  monoWhist_bDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "btagDown",  false, true);
	makehist4(tree, monoWhist_metJetUp.at(itree),  monoWhist_metJetUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jesUp",  false, true);
	makehist4(tree, monoWhist_metJetDw.at(itree),  monoWhist_metJetDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jesDw",  false, true);
	makehist4(tree, monoWhist_metResUp.at(itree),  monoWhist_metResUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jerUp",  false, true);
	makehist4(tree, monoWhist_metResDw.at(itree),  monoWhist_metResDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jerDw",  false, true);
	makehist4(tree, monoWhist_metUncUp.at(itree),  monoWhist_metUncUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "uncUp",  false, true);
	makehist4(tree, monoWhist_metUncDw.at(itree),  monoWhist_metUncDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "uncDw",  false, true);
      }
    }
    itree++;
  }

  itree = 0;
  for(auto tree : monoZtree){
    if(tree){
      cout<<"signal region analysis --> Signal monoZ "<<endl;
      makehist4(tree, monoZhist.at(itree),  monoZhist_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "", false, true);
      if(doShapeSystematics){
	cout<<"signal region analysis --> do signal monoZ sys "<<endl;
	makehist4(tree, monoZhist_bUp.at(itree),  monoZhist_bUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "btagUp",  false, true);
	makehist4(tree, monoZhist_bDw.at(itree),  monoZhist_bDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "btagDown",  false, true);
	makehist4(tree, monoZhist_metJetUp.at(itree),  monoZhist_metJetUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jesUp",  false, true);
	makehist4(tree, monoZhist_metJetDw.at(itree),  monoZhist_metJetDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jesDw",  false, true);
	makehist4(tree, monoZhist_metResUp.at(itree),  monoZhist_metResUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jerUp",  false, true);
	makehist4(tree, monoZhist_metResDw.at(itree),  monoZhist_metResDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "jerDw",  false, true);
	makehist4(tree, monoZhist_metUncUp.at(itree),  monoZhist_metUncUp_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "uncUp",  false, true);
	makehist4(tree, monoZhist_metUncDw.at(itree),  monoZhist_metUncDw_2D.at(itree),  true, 0, category, isWJet, 1.00, lumi, 0, ehists, "uncDw",  false, true);
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
  for(auto sample : monoJhist){ for(auto hist : sample) hist->Write();}
  for(auto sample : monoWhist){ for(auto hist : sample) hist->Write();}  
  for(auto sample : monoZhist){ for(auto hist : sample) hist->Write();}

  if(doShapeSystematics){
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

  for(auto sample : monoJhist_2D){ for (auto hist_2D : sample) hist_2D->Write();}
  for(auto sample : monoWhist_2D){ for (auto hist_2D : sample) hist_2D->Write();}
  for(auto sample : monoZhist_2D){ for (auto hist_2D : sample) hist_2D->Write();}  

  for(auto sample : monoJhist_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
  for(auto sample : monoWhist_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
  for(auto sample : monoZhist_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}

  if(doShapeSystematics){

    for(auto sample : monoJhist_bUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_bDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metJetUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metJetDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metResUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metResDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metUncUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoJhist_metUncDw_2D){ for(auto hist : sample) hist->Write();}

    for(auto sample : monoJhist_bUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_bDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metJetUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metJetDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metResUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metResDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metUncUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoJhist_metUncDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}


    for(auto sample : monoWhist_bUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_bDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metJetUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metJetDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metResUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metResDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metUncUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoWhist_metUncDw_2D){ for(auto hist : sample) hist->Write();}

    for(auto sample : monoWhist_bUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_bDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metJetUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metJetDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metResUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metResDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metUncUp_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}
    for(auto sample : monoWhist_metUncDw_2D){ for (auto hist_2D : sample) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }}

    for(auto sample : monoZhist_bUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_bDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metJetUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metJetDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metResUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metResDw_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metUncUp_2D){ for(auto hist : sample) hist->Write();}
    for(auto sample : monoZhist_metUncDw_2D){ for(auto hist : sample) hist->Write();}

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

void sigdatamchist(TFile* outfile,
                   int category,
                   vector<string> observables,
		   vector<string> observables_2D,
                   double lumi              = 2.24,
                   bool applyQGLReweight    = false,
		   bool doShapeSystematics  = false,
		   bool doAlternativeTop    = false,
                   bool blind = false,
		   bool isHiggsInvisible = false) {

  // Files for Znunu, Wlnu, Zll, top, qcd , diboson, signal, data                                                                                                            
  TFile* znfile  = TFile::Open((baseInputTreePath+"ZJets/sigfilter/sig_tree_ZJetsToNuNu.root").c_str());
  TFile* wlfile  = TFile::Open((baseInputTreePath+"WJets/sigfilter/sig_tree_WJetsToLNu.root").c_str());
  TFile* zlfile  = TFile::Open((baseInputTreePath+"DYJets/sigfilter/sig_tree_DYJetsToLL_M-50.root").c_str());
  TFile* ttfile  = TFile::Open((baseInputTreePath+"Top/sigfilter/sig_tree_Top_amc.root").c_str());
  TFile* qcdfile = TFile::Open((baseInputTreePath+"QCD/sigfilter/sig_tree_QCD.root").c_str());
  TFile* dbfile  = TFile::Open((baseInputTreePath+"DiBoson/sigfilter/sig_tree_DiBoson.root").c_str());
  TFile* gmfile  = TFile::Open((baseInputTreePath+"PhotonJets/sigfilter/sig_tree_GJets.root").c_str());

  // additional top sample
  TFile* ttfile_alt  = NULL;
  if(doAlternativeTop)
    ttfile_alt = TFile::Open((baseInputTreePath+"Top/sigfilter/sig_tree_Top.root").c_str());

  //data                                                                                                                                                                        
  TFile* dtfile = TFile::Open((baseInputTreePath+"MET/sigfilter/sig_tree_crab_MET-Run2015.root").c_str());

  // make met histograms                                                                                                                                                        
  vector<TH1*> znhist;
  vector<TH1*> wlhist;

  vector<TH1*> zlhist;
  vector<TH1*> zlhist_bUp;
  vector<TH1*> zlhist_bDw;
  vector<TH1*> zlhist_metJetUp;
  vector<TH1*> zlhist_metJetDw;
  vector<TH1*> zlhist_metResUp;
  vector<TH1*> zlhist_metResDw;
  vector<TH1*> zlhist_metUncUp;
  vector<TH1*> zlhist_metUncDw;

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

  vector<TH1*> gmhist;
  vector<TH1*> gmhist_bUp;
  vector<TH1*> gmhist_bDw;
  vector<TH1*> gmhist_metJetUp;
  vector<TH1*> gmhist_metJetDw;
  vector<TH1*> gmhist_metResUp;
  vector<TH1*> gmhist_metResDw;
  vector<TH1*> gmhist_metUncUp;
  vector<TH1*> gmhist_metUncDw;

  vector<TH1*> qcdhist;
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
    TH1F* gmhist_temp = new TH1F(("gbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qcdhist_temp = new TH1F(("qbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dthist_temp = new TH1F(("datahist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    znhist.push_back(dynamic_cast<TH1*>(znhist_temp));
    wlhist.push_back(dynamic_cast<TH1*>(wlhist_temp));
    zlhist.push_back(dynamic_cast<TH1*>(zlhist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
    dihist.push_back(dynamic_cast<TH1*>(dihist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
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
      TH1F* gmhist_bUp_temp = new TH1F(("gbkghist_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_bDw_temp = new TH1F(("gbkghist_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_bUp_temp = new TH1F(("zjethist_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_bDw_temp = new TH1F(("zjethist_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      tthist_bUp.push_back(dynamic_cast<TH1*>(tthist_bUp_temp));
      tthist_bDw.push_back(dynamic_cast<TH1*>(tthist_bDw_temp));
      dihist_bUp.push_back(dynamic_cast<TH1*>(dihist_bUp_temp));
      dihist_bDw.push_back(dynamic_cast<TH1*>(dihist_bDw_temp));
      gmhist_bUp.push_back(dynamic_cast<TH1*>(gmhist_bUp_temp));
      gmhist_bDw.push_back(dynamic_cast<TH1*>(gmhist_bDw_temp));
      zlhist_bUp.push_back(dynamic_cast<TH1*>(zlhist_bUp_temp));
      zlhist_bDw.push_back(dynamic_cast<TH1*>(zlhist_bDw_temp));

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


      TH1F* gmhist_metJetUp_temp = new TH1F(("gbkghist_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metJetDw_temp = new TH1F(("gbkghist_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metResUp_temp = new TH1F(("gbkghist_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metResDw_temp = new TH1F(("gbkghist_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metUncUp_temp = new TH1F(("gbkghist_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metUncDw_temp = new TH1F(("gbkghist_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      gmhist_metJetUp.push_back(dynamic_cast<TH1*>(gmhist_metJetUp_temp));
      gmhist_metJetDw.push_back(dynamic_cast<TH1*>(gmhist_metJetDw_temp));
      gmhist_metResUp.push_back(dynamic_cast<TH1*>(gmhist_metResUp_temp));
      gmhist_metResDw.push_back(dynamic_cast<TH1*>(gmhist_metResDw_temp));
      gmhist_metUncUp.push_back(dynamic_cast<TH1*>(gmhist_metUncUp_temp));
      gmhist_metUncDw.push_back(dynamic_cast<TH1*>(gmhist_metUncDw_temp));

      TH1F* zlhist_metJetUp_temp = new TH1F(("zjethist_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_metJetDw_temp = new TH1F(("zjethist_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_metResUp_temp = new TH1F(("zjethist_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_metResDw_temp = new TH1F(("zjethist_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_metUncUp_temp = new TH1F(("zjethist_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_metUncDw_temp = new TH1F(("zjethist_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      zlhist_metJetUp.push_back(dynamic_cast<TH1*>(zlhist_metJetUp_temp));
      zlhist_metJetDw.push_back(dynamic_cast<TH1*>(zlhist_metJetDw_temp));
      zlhist_metResUp.push_back(dynamic_cast<TH1*>(zlhist_metResUp_temp));
      zlhist_metResDw.push_back(dynamic_cast<TH1*>(zlhist_metResDw_temp));
      zlhist_metUncUp.push_back(dynamic_cast<TH1*>(zlhist_metUncUp_temp));
      zlhist_metUncDw.push_back(dynamic_cast<TH1*>(zlhist_metUncDw_temp));
      
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

  vector<TH2*> znhist_2D;
  vector<TH2*> wlhist_2D;
  vector<TH2*> zlhist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_alt_2D;
  vector<TH2*> dihist_2D;
  vector<TH2*> gmhist_2D;
  vector<TH2*> qcdhist_2D;

  vector<TH2*> zlhist_bUp_2D;
  vector<TH2*> zlhist_bDw_2D;
  vector<TH2*> zlhist_metJetUp_2D;
  vector<TH2*> zlhist_metJetDw_2D;
  vector<TH2*> zlhist_metResUp_2D;
  vector<TH2*> zlhist_metResDw_2D;
  vector<TH2*> zlhist_metUncUp_2D;
  vector<TH2*> zlhist_metUncDw_2D;
  vector<TH2*> tthist_bUp_2D;
  vector<TH2*> tthist_bDw_2D;
  vector<TH2*> tthist_metJetUp_2D;
  vector<TH2*> tthist_metJetDw_2D;
  vector<TH2*> tthist_metResUp_2D;
  vector<TH2*> tthist_metResDw_2D;
  vector<TH2*> tthist_metUncUp_2D;
  vector<TH2*> tthist_metUncDw_2D;
  vector<TH2*> tthist_alt_bUp_2D;
  vector<TH2*> tthist_alt_bDw_2D;
  vector<TH2*> tthist_alt_metJetUp_2D;
  vector<TH2*> tthist_alt_metJetDw_2D;
  vector<TH2*> tthist_alt_metResUp_2D;
  vector<TH2*> tthist_alt_metResDw_2D;
  vector<TH2*> tthist_alt_metUncUp_2D;
  vector<TH2*> tthist_alt_metUncDw_2D;
  vector<TH2*> dihist_bUp_2D;
  vector<TH2*> dihist_bDw_2D;
  vector<TH2*> dihist_metJetUp_2D;
  vector<TH2*> dihist_metJetDw_2D;
  vector<TH2*> dihist_metResUp_2D;
  vector<TH2*> dihist_metResDw_2D;
  vector<TH2*> dihist_metUncUp_2D;
  vector<TH2*> dihist_metUncDw_2D;
  vector<TH2*> gmhist_bUp_2D;
  vector<TH2*> gmhist_bDw_2D;
  vector<TH2*> gmhist_metJetUp_2D;
  vector<TH2*> gmhist_metJetDw_2D;
  vector<TH2*> gmhist_metResUp_2D;
  vector<TH2*> gmhist_metResDw_2D;
  vector<TH2*> gmhist_metUncUp_2D;
  vector<TH2*> gmhist_metUncDw_2D;
  vector<TH2*> dthist_2D;


  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
    
    TH2F* znhist_temp = new TH2F(("zinvhist_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* wlhist_temp = new TH2F(("wjethist_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* zlhist_temp = new TH2F(("zjethist_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* tthist_temp = new TH2F(("tbkghist_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* dihist_temp = new TH2F(("dbkghist_2D"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* gmhist_temp = new TH2F(("gbkghist_2D"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* qcdhist_temp = new TH2F(("qbkghist_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* dthist_temp = new TH2F(("datahist_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

    znhist_2D.push_back(dynamic_cast<TH2*>(znhist_temp));
    wlhist_2D.push_back(dynamic_cast<TH2*>(wlhist_temp));
    zlhist_2D.push_back(dynamic_cast<TH2*>(zlhist_temp));
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    qcdhist_2D.push_back(dynamic_cast<TH2*>(qcdhist_temp));
    dihist_2D.push_back(dynamic_cast<TH2*>(dihist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));

    if(doAlternativeTop){
      TH2F* tthist_alt_temp = new TH2F(("tbkghist_alt_2D_"+obs).c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      tthist_alt_2D.push_back(dynamic_cast<TH2*>(tthist_alt_temp));
    }    
  }

  
  TTree* zntree = (TTree*)znfile->Get("tree/tree");
  TTree* wltree = (TTree*)wlfile->Get("tree/tree");
  TTree* zltree = (TTree*)zlfile->Get("tree/tree");
  TTree* tttree = (TTree*)ttfile->Get("tree/tree");
  TTree* tttree_alt = NULL;
  if(ttfile_alt and doAlternativeTop)
    tttree_alt = (TTree*)ttfile_alt->Get("tree/tree");

  TTree* ditree  = (TTree*)dbfile->Get("tree/tree");
  TTree* gmtree  = (TTree*)gmfile->Get("tree/tree");
  TTree* qcdtree = (TTree*)qcdfile->Get("tree/tree");

  TTree* dttree = (TTree*)dtfile->Get("tree");

  // get k-factors NLO                                                                                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  if(not znlohist or not zlohist or not zewkhist){
    znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    zlohist  = (TH1*)kffile.Get("zlo/zlo_nominal");
    zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
    znlohist->Divide(zlohist);
  }

  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");
 
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  if(not wnlohist or not wlohist or not wewkhist){
    wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    wlohist  = (TH1*)kffile.Get("wlo/wlo_nominal");
    wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");
    wnlohist->Divide(wlohist);
  }

  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  if(not anlohist or not alohist or not aewkhist){
    anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
    alohist  = (TH1*)kffile.Get("alo/alo_nominal");
    aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");
    anlohist->Divide(alohist);
  }

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> whists;
  vector<TH1*> ahists;

  // apply EWK and QCD corrections                                                                                                                                              
  zhists.push_back(znlohist); 
  zhists.push_back(zewkhist);
  whists.push_back(wnlohist); 
  whists.push_back(wewkhist);    
  ahists.push_back(anlohist);
  ahists.push_back(aewkhist);


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
  makehist4(zntree,znhist,znhist_2D,true,0,category,false,1.00,lumi,QGLZ_index,zhists,"",false,true,0,isHiggsInvisible);
  cout<<"signal region: W+jets sample "<<endl;
  makehist4(wltree,wlhist,wlhist_2D,true,0,category,false,1.00,lumi,QGLW_index,whists,"",false,true,0,isHiggsInvisible);
  cout<<"signal region: Z+jets sample "<<endl;
  makehist4(zltree,zlhist,zlhist_2D,true,0,category,false,1.00,lumi,QGLZ_index,zhists,"",false,true,0,isHiggsInvisible);
  cout<<"signal region: gamma+jets sample "<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,0,category,false,1.00,lumi,QGLG_index,ahists,"",false,true,0,isHiggsInvisible);
  cout<<"signal region: TTbar sample "<<endl;
  makehist4(tttree,tthist,tthist_2D,true,0,category,false,1.00,lumi,QGLT_index,ehists,"",true,true,0,isHiggsInvisible);

    //alternative ttbar             
  if(doAlternativeTop){
    cout<<"signal region: TTbar alternative sample "<<endl;
    makehist4(tttree_alt, tthist_alt,  tthist_alt_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "", false, true,0,isHiggsInvisible);
  }

  cout<<"signal region: Diboson sample "<<endl;
  makehist4(ditree, dihist,  dihist_2D,  true, 0, category, isWJet, 1.00, lumi,   0, ehists, "",  false, true);
  cout<<"signal region: QCD sample "<<endl;
  makehist4(qcdtree, qcdhist,  qcdhist_2D,  true, 0, category, false, 1.00, lumi, 0,  ehists, "", false, true);

  if(doShapeSystematics){
    cout<<"signal region analysis --> do top shape sys "<<endl;
    makehist4(tttree, tthist_bUp,  tthist_bUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "btagUp",   true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_bDw,  tthist_bDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "btagDown", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metJetUp,  tthist_metJetUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jesUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metJetDw,  tthist_metJetDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jesDw", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metResUp,  tthist_metResUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jerUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metResDw,  tthist_metResDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "jerDw", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metUncUp,  tthist_metUncUp_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "uncUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metUncDw,  tthist_metUncDw_2D,  true, 0, category, false, 1.00, lumi,    QGLT_index, ehists, "uncDw", true, true,0,isHiggsInvisible);
    
    if(doAlternativeTop){
      cout<<"signal region analysis --> do top alternative shape sys "<<endl;
      makehist4(tttree_alt, tthist_alt_bUp, tthist_alt_bUp_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "btagUp", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_bDw, tthist_alt_bDw_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "btagDown", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_metJetUp,  tthist_alt_metJetUp_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "jesUp", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_metJetDw,  tthist_alt_metJetDw_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "jesDw", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_metResUp,  tthist_alt_metResUp_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "jerUp", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_metResDw,  tthist_alt_metResDw_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "jerDw", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_metUncUp,  tthist_alt_metUncUp_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "uncUp", true, true,0,isHiggsInvisible);
      makehist4(tttree_alt, tthist_alt_metUncDw,  tthist_alt_metUncDw_2D,  true, 0, category, false, 1.00, lumi,QGLT_index, ehists, "uncDw", true, true,0,isHiggsInvisible);
    }


    cout<<"signal region analysis --> do diboson shape sys "<<endl;
    makehist4(ditree, dihist_bUp,  dihist_bUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "btagUp", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_bDw,  dihist_bDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "btagDown", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_metJetUp,  dihist_metJetUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jesUp", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_metJetDw,  dihist_metJetDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jesDw", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_metResUp,  dihist_metResUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jerUp", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_metResDw,  dihist_metResDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "jerDw", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_metUncUp,  dihist_metUncUp_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "uncUp", false, true, 0,isHiggsInvisible);
    makehist4(ditree, dihist_metUncDw,  dihist_metUncDw_2D,  true, 0, category, isWJet, 1.00, lumi,    0, ehists, "uncDw", false, true, 0,isHiggsInvisible);

    cout<<"signal region analysis --> do DYJets shape sys "<<endl;
    makehist4(zltree, zlhist_bUp,  zlhist_bUp_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "btagUp", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_bDw,  zlhist_bDw_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "btagDown", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_metJetUp,  zlhist_metJetUp_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "jesUp", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_metJetDw,  zlhist_metJetDw_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "jesDw", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_metResUp,  zlhist_metResUp_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "jerUp", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_metResDw,  zlhist_metResDw_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "jerDw", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_metUncUp,  zlhist_metUncUp_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "uncUp", false, true, 0,isHiggsInvisible);
    makehist4(zltree, zlhist_metUncDw,  zlhist_metUncDw_2D,  true, 0, category, false, 1.00, lumi,    QGLZ_index, zhists, "uncDw", false, true, 0,isHiggsInvisible);

    cout<<"signal region analysis --> do gamma+jets shape sys "<<endl;
    makehist4(gmtree, gmhist_bUp,  gmhist_bUp_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "btagUp", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_bDw,  gmhist_bDw_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "btagDown", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metJetUp,  gmhist_metJetUp_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "jesUp", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metJetDw,  gmhist_metJetDw_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "jesDw", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metResUp,  gmhist_metResUp_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "jerUp", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metResDw,  gmhist_metResDw_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "jerDw", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metUncUp,  gmhist_metUncUp_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "uncUp", false, true, 0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metUncDw,  gmhist_metUncDw_2D,  true, 0, category, false, 1.00, lumi,    QGLG_index, ahists, "uncDw", false, true, 0,isHiggsInvisible);
    
  }

  // take average of ttbar                                                                                                                                                     
  if(doAlternativeTop){

    for(size_t iHisto = 0; iHisto < tthist.size(); iHisto++)
      makeAverage(tthist.at(iHisto),tthist_alt.at(iHisto));      

    for(size_t iHisto = 0; iHisto < tthist_2D.size(); iHisto++)
      makeAverage(tthist_2D.at(iHisto),tthist_alt_2D.at(iHisto));

    if(doShapeSystematics)
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

  //smooth
  for(auto hist : tthist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dihist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : gmhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dihist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : zlhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : qcdhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  for(auto hist : tthist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dihist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : gmhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dihist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : zlhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : qcdhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  if(doShapeSystematics){
    
    for(auto hist : tthist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : dihist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : zlhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    
    for(auto hist : tthist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : dihist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : zlhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
            
  }

  // data                                                
  cout<<"signal region analysis --> loop on data "<<endl;
  makehist4(dttree, dthist,  dthist_2D, false, 0, category, false, 1.00, lumi, 0, ehists, "", false, true, 0,isHiggsInvisible);
    
  if (blind) {
    for( size_t ihist = 0; ihist < dthist.size(); ihist++){
      for (int i = 1; i <= dthist.at(ihist)->GetNbinsX(); i++) {
        double binval = 0.0;
        binval += znhist.at(ihist)->GetBinContent(i);
	binval += wlhist.at(ihist)->GetBinContent(i);
        binval += zlhist.at(ihist)->GetBinContent(i);
        binval += tthist.at(ihist)->GetBinContent(i);
        binval += dihist.at(ihist)->GetBinContent(i);
        binval += gmhist.at(ihist)->GetBinContent(i);
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
          binval += gmhist_2D.at(ihist)->GetBinContent(iX,iY);
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
  for(auto hist : gmhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();
  for(auto hist : dthist) hist->Write();

  //
  if(doShapeSystematics){
    
    for(auto hist : zlhist_bUp) hist->Write();
    for(auto hist : zlhist_bDw) hist->Write();
    for(auto hist : zlhist_metJetUp) hist->Write();
    for(auto hist : zlhist_metJetDw) hist->Write();
    for(auto hist : zlhist_metResUp) hist->Write();
    for(auto hist : zlhist_metResDw) hist->Write();
    for(auto hist : zlhist_metUncUp) hist->Write();
    for(auto hist : zlhist_metUncDw) hist->Write();

    for(auto hist : gmhist_bUp) hist->Write();
    for(auto hist : gmhist_bDw) hist->Write();
    for(auto hist : gmhist_metJetUp) hist->Write();
    for(auto hist : gmhist_metJetDw) hist->Write();
    for(auto hist : gmhist_metResUp) hist->Write();
    for(auto hist : gmhist_metResDw) hist->Write();
    for(auto hist : gmhist_metUncUp) hist->Write();
    for(auto hist : gmhist_metUncDw) hist->Write();

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
    
  }

  // store hist_2Dograms                                                                                                                                                     
  for(auto hist_2D : znhist_2D) hist_2D->Write();
  for(auto hist_2D : wlhist_2D) hist_2D->Write();
  for(auto hist_2D : zlhist_2D) hist_2D->Write();
  for(auto hist_2D : gmhist_2D) hist_2D->Write();
  for(auto hist_2D : tthist_2D) hist_2D->Write();
  for(auto hist_2D : dihist_2D) hist_2D->Write();
  for(auto hist_2D : qcdhist_2D) hist_2D->Write();
  for(auto hist_2D : dthist_2D) hist_2D->Write();

  for(auto hist_2D : znhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : wlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : zlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : gmhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : tthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dihist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : qcdhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }


  if(doShapeSystematics){
    
    for(auto hist : zlhist_bUp_2D) hist->Write();
    for(auto hist : zlhist_bDw_2D) hist->Write();
    for(auto hist : zlhist_metJetUp_2D) hist->Write();
    for(auto hist : zlhist_metJetDw_2D) hist->Write();
    for(auto hist : zlhist_metResUp_2D) hist->Write();
    for(auto hist : zlhist_metResDw_2D) hist->Write();
    for(auto hist : zlhist_metUncUp_2D) hist->Write();
    for(auto hist : zlhist_metUncDw_2D) hist->Write();

    for(auto hist_2D : zlhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    for(auto hist : gmhist_bUp_2D) hist->Write();
    for(auto hist : gmhist_bDw_2D) hist->Write();
    for(auto hist : gmhist_metJetUp_2D) hist->Write();
    for(auto hist : gmhist_metJetDw_2D) hist->Write();
    for(auto hist : gmhist_metResUp_2D) hist->Write();
    for(auto hist : gmhist_metResDw_2D) hist->Write();
    for(auto hist : gmhist_metUncUp_2D) hist->Write();
    for(auto hist : gmhist_metUncDw_2D) hist->Write();

    for(auto hist_2D : gmhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    for(auto hist : tthist_bUp_2D) hist->Write();
    for(auto hist : tthist_bDw_2D) hist->Write();
    for(auto hist : tthist_metJetUp_2D) hist->Write();
    for(auto hist : tthist_metJetDw_2D) hist->Write();
    for(auto hist : tthist_metResUp_2D) hist->Write();
    for(auto hist : tthist_metResDw_2D) hist->Write();
    for(auto hist : tthist_metUncUp_2D) hist->Write();
    for(auto hist : tthist_metUncDw_2D) hist->Write();

    for(auto hist_2D : tthist_bUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_bDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    for(auto hist : dihist_bUp_2D) hist->Write();
    for(auto hist : dihist_bDw_2D) hist->Write();
    for(auto hist : dihist_metJetUp_2D) hist->Write();
    for(auto hist : dihist_metJetDw_2D) hist->Write();
    for(auto hist : dihist_metResUp_2D) hist->Write();
    for(auto hist : dihist_metResDw_2D) hist->Write();
    for(auto hist : dihist_metUncUp_2D) hist->Write();
    for(auto hist : dihist_metUncDw_2D) hist->Write();

    for(auto hist_2D : dihist_bUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_bDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    
  }

  znfile->Close();
  wlfile->Close();
  zlfile->Close();
  ttfile->Close();
  dbfile->Close();
  gmfile->Close();
  qcdfile->Close();
  dtfile->Close();
  kffile.Close();
  if(ttfile_alt) ttfile_alt->Close();

  znhist.clear();
  wlhist.clear();
  zlhist.clear();
  zlhist_bUp.clear();
  zlhist_bDw.clear();
  zlhist_metJetUp.clear();
  zlhist_metJetDw.clear();
  zlhist_metResUp.clear();
  zlhist_metResDw.clear();
  zlhist_metUncUp.clear();
  zlhist_metUncDw.clear();
  tthist.clear();
  tthist_bUp.clear();
  tthist_bDw.clear();
  tthist_metJetUp.clear();
  tthist_metJetDw.clear();
  tthist_metResUp.clear();
  tthist_metResDw.clear();
  tthist_metUncUp.clear();
  tthist_metUncDw.clear();
  tthist_alt.clear();
  tthist_alt_bUp.clear();
  tthist_alt_bDw.clear();
  tthist_alt_metJetUp.clear();
  tthist_alt_metJetDw.clear();
  tthist_alt_metResUp.clear();
  tthist_alt_metResDw.clear();
  tthist_alt_metUncUp.clear();
  tthist_alt_metUncDw.clear();
  dihist.clear();
  dihist_bUp.clear();
  dihist_bDw.clear();
  dihist_metJetUp.clear();
  dihist_metJetDw.clear();
  dihist_metResUp.clear();
  dihist_metResDw.clear();
  dihist_metUncUp.clear();
  dihist_metUncDw.clear();
  gmhist.clear();
  gmhist_bUp.clear();
  gmhist_bDw.clear();
  gmhist_metJetUp.clear();
  gmhist_metJetDw.clear();
  gmhist_metResUp.clear();
  gmhist_metResDw.clear();
  gmhist_metUncUp.clear();
  gmhist_metUncDw.clear();
  qcdhist.clear();
  dthist.clear();

  znhist_2D.clear();
  wlhist_2D.clear();
  zlhist_2D.clear();
  zlhist_bUp_2D.clear();
  zlhist_bDw_2D.clear();
  zlhist_metJetUp_2D.clear();
  zlhist_metJetDw_2D.clear();
  zlhist_metResUp_2D.clear();
  zlhist_metResDw_2D.clear();
  zlhist_metUncUp_2D.clear();
  zlhist_metUncDw_2D.clear();
  tthist_2D.clear();
  tthist_bUp_2D.clear();
  tthist_bDw_2D.clear();
  tthist_metJetUp_2D.clear();
  tthist_metJetDw_2D.clear();
  tthist_metResUp_2D.clear();
  tthist_metResDw_2D.clear();
  tthist_metUncUp_2D.clear();
  tthist_metUncDw_2D.clear();
  tthist_alt_2D.clear();
  tthist_alt_bUp_2D.clear();
  tthist_alt_bDw_2D.clear();
  tthist_alt_metJetUp_2D.clear();
  tthist_alt_metJetDw_2D.clear();
  tthist_alt_metResUp_2D.clear();
  tthist_alt_metResDw_2D.clear();
  tthist_alt_metUncUp_2D.clear();
  tthist_alt_metUncDw_2D.clear();
  dihist_2D.clear();
  dihist_bUp_2D.clear();
  dihist_bDw_2D.clear();
  dihist_metJetUp_2D.clear();
  dihist_metJetDw_2D.clear();
  dihist_metResUp_2D.clear();
  dihist_metResDw_2D.clear();
  dihist_metUncUp_2D.clear();
  dihist_metUncDw_2D.clear();
  gmhist_2D.clear();
  gmhist_bUp_2D.clear();
  gmhist_bDw_2D.clear();
  gmhist_metJetUp_2D.clear();
  gmhist_metJetDw_2D.clear();
  gmhist_metResUp_2D.clear();
  gmhist_metResDw_2D.clear();
  gmhist_metUncUp_2D.clear();
  gmhist_metUncDw_2D.clear();
  qcdhist_2D.clear();
  dthist_2D.clear();


  cout << "Templates for the signal region computed ..." << endl;
}


// build templates for photon+jets control region                                                                                                                           
void gamdatamchist(TFile* outfile,
                  int category,
                   vector<string> observables,
                   vector<string> observables_2D,
                   double lumi           = 2.24,
                   bool applyQGLReweight = false,
		   bool isHiggsInvisible = false
                   ) {


  TFile* dtfile =  TFile::Open((baseInputTreePath+"SinglePhoton/gamfilter/gam_tree_crab_SinglePhoton-Run2015.root").c_str());
  TFile* gmfile =  TFile::Open((baseInputTreePath+"PhotonJets/gamfilter/gam_tree_GJets.root").c_str());

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

 
  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty() )
      cout<<"No binning for this observable --> please define it"<<endl;

    TH2F* gmhist_temp = new TH2F(("gbkghistgam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* qchist_temp = new TH2F(("qbkghistgam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* dthist_temp = new TH2F(("datahistgam_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

    qcdhist_2D.push_back(dynamic_cast<TH2*>(qchist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));
  }


  TTree* dttree = (TTree*)dtfile->Get("tree");
  TTree* gmtree = (TTree*)gmfile->Get("tree/tree");

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                           
  // get k-factors NLO                                                                                                                                                        
  TFile kffile (kfactorFile.c_str());
  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);
  
  if(not anlohist or not alohist or not aewkhist){
    anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
    alohist  = (TH1*)kffile.Get("alo/alo_nominal");
    aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");
    anlohist->Divide(alohist);
  }

  vector<TH1*> ahists;
  vector<TH1*> ehists;

  ahists.push_back(anlohist);
  ahists.push_back(aewkhist);

  if(applyQGLReweight){
    cout<<"gamma+jets control region --> data"<<endl;
    makehist4(dttree, dthist, dthist_2D, false, 5, category, false, 1.00, lumi, 0, ehists, "", false, true, 0,isHiggsInvisible);
    cout<<"gamma+jets control region --> gamma+jets"<<endl;
    makehist4(gmtree, gmhist, gmhist_2D, true,  5, category, false, 1.00, lumi, 3, ahists, "", false, true, 0,isHiggsInvisible);
    cout<<"gamma+jets control region: QCD"<<endl;
    makehist4(dttree, qcdhist, qcdhist_2D, false, 6, category, false, 1.00, lumi, 0, ehists, "", false, true, 0,isHiggsInvisible);
  }
  else{
    cout<<"gamma+jets control region --> data"<<endl;
    makehist4(dttree, dthist, dthist_2D, false, 5, category, false, 1.00, lumi, 0, ehists, "", false, true, 0,isHiggsInvisible);
    cout<<"gamma+jets control region --> gamma+jets"<<endl;
    makehist4(gmtree, gmhist, gmhist_2D, true,  5, category, false, 1.00, lumi, 0, ahists, "", false, true, 0,isHiggsInvisible);
    cout<<"gamma+jets control region --> QCD"<<endl;
    makehist4(dttree, qcdhist, qcdhist_2D, false, 6, category, false, 1.00, lumi, 0, ehists, "", false, true, 0,isHiggsInvisible);
  }
			    
  outfile->cd();

  for(auto hist : dthist) hist->Write();
  for(auto hist : gmhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();

  for(auto hist : dthist_2D) hist->Write();
  for(auto hist : gmhist_2D) hist->Write();
  for(auto hist : qcdhist_2D) hist->Write();

  for(auto hist : dthist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
  for(auto hist : gmhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
  for(auto hist : qcdhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

  dtfile->Close();
  gmfile->Close();
  kffile.Close();

  dthist.clear();
  qcdhist.clear();
  gmhist.clear();
  dthist_2D.clear();
  qcdhist_2D.clear();
  gmhist_2D.clear();

  cout << "Templates for the gamma+jets control region computed ..." << endl;
}


//build templates for Zmumu, Zee, Wenu, Wmunu                                                                                                                                  
void lepdatamchist(TFile* outfile, int sample, int category, vector<string> observables, vector<string> observables_2D,
		   double lumi = 2.24, bool applyQGLReweight = false, bool doShapeSystematics = false, bool isHiggsInvisible = false) {

  if (sample != 1 && sample != 2 && sample != 3 && sample != 4) return;

  TFile* ttfile  = NULL;
  TFile* dbfile  = NULL; 
  TFile* gmfile  = NULL;
  TFile* qcfile  = NULL;
  TFile* vlfile  = NULL;
  TFile* vllfile = NULL;
  TFile* dtfile  = NULL;
  TFile* dtfile_2  = NULL;

  string suffix;

  if(sample == 1){

    suffix = "zmm";
    vllfile = TFile::Open((baseInputTreePath+"DYJets/zmmfilter/zmm_tree_DYJetsToLL_M-50.root").c_str());
    vlfile = TFile::Open((baseInputTreePath+"WJets/zmmfilter/zmm_tree_WJetsToLNu.root").c_str());
    qcfile = TFile::Open((baseInputTreePath+"QCD/zmmfilter/zmm_tree_QCD.root").c_str());
    dbfile = TFile::Open((baseInputTreePath+"DiBoson/zmmfilter/zmm_tree_DiBoson.root").c_str());
    gmfile = TFile::Open((baseInputTreePath+"PhotonJets/zmmfilter/zmm_tree_GJets.root").c_str());
    ttfile = TFile::Open((baseInputTreePath+"Top/zmmfilter/zmm_tree_Top_amc.root").c_str());
    dtfile = TFile::Open((baseInputTreePath+"MET/zmmfilter/zmm_tree_crab_MET-Run2015.root").c_str());
  }
  else if(sample == 2){

    suffix = "wmn";
    vllfile = TFile::Open((baseInputTreePath+"DYJets/wmnfilter/wmn_tree_DYJetsToLL_M-50.root").c_str());
    vlfile = TFile::Open((baseInputTreePath+"WJets/wmnfilter/wmn_tree_WJetsToLNu.root").c_str());
    qcfile = TFile::Open((baseInputTreePath+"QCD/wmnfilter/wmn_tree_QCD.root").c_str());
    dbfile = TFile::Open((baseInputTreePath+"DiBoson/wmnfilter/wmn_tree_DiBoson.root").c_str());
    gmfile = TFile::Open((baseInputTreePath+"PhotonJets/wmnfilter/wmn_tree_GJets.root").c_str());
    ttfile = TFile::Open((baseInputTreePath+"Top/wmnfilter/wmn_tree_Top_amc.root").c_str());
    dtfile = TFile::Open((baseInputTreePath+"MET/wmnfilter/wmn_tree_crab_MET-Run2015.root").c_str());
  }
  else if(sample == 3){

    suffix = "zee";
    vllfile = TFile::Open((baseInputTreePath+"DYJets/zeefilter/zee_tree_DYJetsToLL_M-50.root").c_str());
    vlfile = TFile::Open((baseInputTreePath+"WJets/zeefilter/zee_tree_WJetsToLNu.root").c_str());
    qcfile = TFile::Open((baseInputTreePath+"QCD/zeefilter/zee_tree_QCD.root").c_str());
    dbfile = TFile::Open((baseInputTreePath+"DiBoson/zeefilter/zee_tree_DiBoson.root").c_str());
    gmfile = TFile::Open((baseInputTreePath+"PhotonJets/zeefilter/zee_tree_GJets.root").c_str());
    ttfile = TFile::Open((baseInputTreePath+"Top/zeefilter/zee_tree_Top_amc.root").c_str());
    dtfile = TFile::Open((baseInputTreePath+"SingleElectron/zeefilter/zee_tree_crab_SingleEle-Run2015.root").c_str());
    dtfile_2 = TFile::Open((baseInputTreePath+"SinglePhoton/zeefilter/zee_tree_crab_SinglePhoton-Run2015.root").c_str());
  }
  else if(sample == 4){

    suffix = "wen";
    vllfile = TFile::Open((baseInputTreePath+"DYJets/wenfilter/wen_tree_DYJetsToLL_M-50.root").c_str());
    vlfile = TFile::Open((baseInputTreePath+"WJets/wenfilter/wen_tree_WJetsToLNu.root").c_str());
    qcfile = TFile::Open((baseInputTreePath+"QCD/wenfilter/wen_tree_QCD.root").c_str());
    dbfile = TFile::Open((baseInputTreePath+"DiBoson/wenfilter/wen_tree_DiBoson.root").c_str());
    gmfile = TFile::Open((baseInputTreePath+"PhotonJets/wenfilter/wen_tree_GJets.root").c_str());
    ttfile = TFile::Open((baseInputTreePath+"Top/wenfilter/wen_tree_Top_amc.root").c_str());
    dtfile = TFile::Open((baseInputTreePath+"SingleElectron/wenfilter/wen_tree_crab_SingleEle-Run2015.root").c_str());
    dtfile_2 = TFile::Open((baseInputTreePath+"SinglePhoton/wenfilter/wen_tree_crab_SinglePhoton-Run2015.root").c_str());
  }


  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_bUp;
  vector<TH1*> tthist_bDw;
  vector<TH1*> tthist_metJetUp;
  vector<TH1*> tthist_metJetDw;
  vector<TH1*> tthist_metResUp;
  vector<TH1*> tthist_metResDw;
  vector<TH1*> tthist_metUncUp;
  vector<TH1*> tthist_metUncDw;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> dbhist_bUp;
  vector<TH1*> dbhist_bDw;
  vector<TH1*> dbhist_metJetUp;
  vector<TH1*> dbhist_metJetDw;
  vector<TH1*> dbhist_metResUp;
  vector<TH1*> dbhist_metResDw;
  vector<TH1*> dbhist_metUncUp;
  vector<TH1*> dbhist_metUncDw;
  vector<TH1*> gmhist;
  vector<TH1*> gmhist_bUp;
  vector<TH1*> gmhist_bDw;
  vector<TH1*> gmhist_metJetUp;
  vector<TH1*> gmhist_metJetDw;
  vector<TH1*> gmhist_metResUp;
  vector<TH1*> gmhist_metResDw;
  vector<TH1*> gmhist_metUncUp;
  vector<TH1*> gmhist_metUncDw;
  vector<TH1*> vlhist;
  vector<TH1*> vlhist_bUp;
  vector<TH1*> vlhist_bDw;
  vector<TH1*> vlhist_metJetUp;
  vector<TH1*> vlhist_metJetDw;
  vector<TH1*> vlhist_metResUp;
  vector<TH1*> vlhist_metResDw;
  vector<TH1*> vlhist_metUncUp;
  vector<TH1*> vlhist_metUncDw;
  vector<TH1*> vllhist;
  vector<TH1*> vllhist_bUp;
  vector<TH1*> vllhist_bDw;
  vector<TH1*> vllhist_metJetUp;
  vector<TH1*> vllhist_metJetDw;
  vector<TH1*> vllhist_metResUp;
  vector<TH1*> vllhist_metResDw;
  vector<TH1*> vllhist_metUncUp;
  vector<TH1*> vllhist_metUncDw;

  vector<float> bins;
  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;


    TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* gmhist_temp = new TH1F((string("gbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));

    if(doShapeSystematics){

      TH1F* tthist_bUp_temp = new TH1F((string("tbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_bDw_temp = new TH1F((string("tbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metJetUp_temp = new TH1F((string("tbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metJetDw_temp = new TH1F((string("tbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metResUp_temp = new TH1F((string("tbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metResDw_temp = new TH1F((string("tbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metUncUp_temp = new TH1F((string("tbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_metUncDw_temp = new TH1F((string("tbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      tthist_bUp.push_back(dynamic_cast<TH1*>(tthist_bUp_temp));
      tthist_bDw.push_back(dynamic_cast<TH1*>(tthist_bDw_temp));
      tthist_metJetUp.push_back(dynamic_cast<TH1*>(tthist_metJetUp_temp));
      tthist_metJetDw.push_back(dynamic_cast<TH1*>(tthist_metJetDw_temp));
      tthist_metResUp.push_back(dynamic_cast<TH1*>(tthist_metResUp_temp));
      tthist_metResDw.push_back(dynamic_cast<TH1*>(tthist_metResDw_temp));
      tthist_metUncUp.push_back(dynamic_cast<TH1*>(tthist_metUncUp_temp));
      tthist_metUncDw.push_back(dynamic_cast<TH1*>(tthist_metUncDw_temp));

      TH1F* dbhist_bUp_temp = new TH1F((string("dbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_bDw_temp = new TH1F((string("dbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metJetUp_temp = new TH1F((string("dbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metJetDw_temp = new TH1F((string("dbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metResUp_temp = new TH1F((string("dbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metResDw_temp = new TH1F((string("dbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metUncUp_temp = new TH1F((string("dbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metUncDw_temp = new TH1F((string("dbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      dbhist_bUp.push_back(dynamic_cast<TH1*>(dbhist_bUp_temp));
      dbhist_bDw.push_back(dynamic_cast<TH1*>(dbhist_bDw_temp));
      dbhist_metJetUp.push_back(dynamic_cast<TH1*>(dbhist_metJetUp_temp));
      dbhist_metJetDw.push_back(dynamic_cast<TH1*>(dbhist_metJetDw_temp));
      dbhist_metResUp.push_back(dynamic_cast<TH1*>(dbhist_metResUp_temp));
      dbhist_metResDw.push_back(dynamic_cast<TH1*>(dbhist_metResDw_temp));
      dbhist_metUncUp.push_back(dynamic_cast<TH1*>(dbhist_metUncUp_temp));
      dbhist_metUncDw.push_back(dynamic_cast<TH1*>(dbhist_metUncDw_temp));

      TH1F* gmhist_bUp_temp = new TH1F((string("gbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_bDw_temp = new TH1F((string("gbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metJetUp_temp = new TH1F((string("gbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metJetDw_temp = new TH1F((string("gbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metResUp_temp = new TH1F((string("gbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metResDw_temp = new TH1F((string("gbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metUncUp_temp = new TH1F((string("gbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metUncDw_temp = new TH1F((string("gbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      gmhist_bUp.push_back(dynamic_cast<TH1*>(gmhist_bUp_temp));
      gmhist_bDw.push_back(dynamic_cast<TH1*>(gmhist_bDw_temp));
      gmhist_metJetUp.push_back(dynamic_cast<TH1*>(gmhist_metJetUp_temp));
      gmhist_metJetDw.push_back(dynamic_cast<TH1*>(gmhist_metJetDw_temp));
      gmhist_metResUp.push_back(dynamic_cast<TH1*>(gmhist_metResUp_temp));
      gmhist_metResDw.push_back(dynamic_cast<TH1*>(gmhist_metResDw_temp));
      gmhist_metUncUp.push_back(dynamic_cast<TH1*>(gmhist_metUncUp_temp));
      gmhist_metUncDw.push_back(dynamic_cast<TH1*>(gmhist_metUncDw_temp));

    }

    if((sample == 1 or sample == 3) and doShapeSystematics){
      TH1F* vlhist_bUp_temp = new TH1F((string("vlbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_bDw_temp = new TH1F((string("vlbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metJetUp_temp = new TH1F((string("vlbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metJetDw_temp = new TH1F((string("vlbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metResUp_temp = new TH1F((string("vlbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metResDw_temp = new TH1F((string("vlbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metUncUp_temp = new TH1F((string("vlbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metUncDw_temp = new TH1F((string("vlbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      vlhist_bUp.push_back(dynamic_cast<TH1*>(vlhist_bUp_temp));
      vlhist_bDw.push_back(dynamic_cast<TH1*>(vlhist_bDw_temp));
      vlhist_metJetUp.push_back(dynamic_cast<TH1*>(vlhist_metJetUp_temp));
      vlhist_metJetDw.push_back(dynamic_cast<TH1*>(vlhist_metJetDw_temp));
      vlhist_metResUp.push_back(dynamic_cast<TH1*>(vlhist_metResUp_temp));
      vlhist_metResDw.push_back(dynamic_cast<TH1*>(vlhist_metResDw_temp));
      vlhist_metUncUp.push_back(dynamic_cast<TH1*>(vlhist_metUncUp_temp));
      vlhist_metUncDw.push_back(dynamic_cast<TH1*>(vlhist_metUncDw_temp));
      
    }
    else if((sample == 2 or sample == 4) and doShapeSystematics){

      TH1F* vllhist_bUp_temp = new TH1F((string("vllbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_bDw_temp = new TH1F((string("vllbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metJetUp_temp = new TH1F((string("vllbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metJetDw_temp = new TH1F((string("vllbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metResUp_temp = new TH1F((string("vllbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metResDw_temp = new TH1F((string("vllbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metUncUp_temp = new TH1F((string("vllbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metUncDw_temp = new TH1F((string("vllbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      vllhist_bUp.push_back(dynamic_cast<TH1*>(vllhist_bUp_temp));
      vllhist_bDw.push_back(dynamic_cast<TH1*>(vllhist_bDw_temp));
      vllhist_metJetUp.push_back(dynamic_cast<TH1*>(vllhist_metJetUp_temp));
      vllhist_metJetDw.push_back(dynamic_cast<TH1*>(vllhist_metJetDw_temp));
      vllhist_metResUp.push_back(dynamic_cast<TH1*>(vllhist_metResUp_temp));
      vllhist_metResDw.push_back(dynamic_cast<TH1*>(vllhist_metResDw_temp));
      vllhist_metUncUp.push_back(dynamic_cast<TH1*>(vllhist_metUncUp_temp));
      vllhist_metUncDw.push_back(dynamic_cast<TH1*>(vllhist_metUncDw_temp));
     
    }

  }

  vector<TH2*> dthist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_bUp_2D;
  vector<TH2*> tthist_bDw_2D;
  vector<TH2*> tthist_metJetUp_2D;
  vector<TH2*> tthist_metJetDw_2D;
  vector<TH2*> tthist_metResUp_2D;
  vector<TH2*> tthist_metResDw_2D;
  vector<TH2*> tthist_metUncUp_2D;
  vector<TH2*> tthist_metUncDw_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> dbhist_bUp_2D;
  vector<TH2*> dbhist_bDw_2D;
  vector<TH2*> dbhist_metJetUp_2D;
  vector<TH2*> dbhist_metJetDw_2D;
  vector<TH2*> dbhist_metResUp_2D;
  vector<TH2*> dbhist_metResDw_2D;
  vector<TH2*> dbhist_metUncUp_2D;
  vector<TH2*> dbhist_metUncDw_2D; 
  vector<TH2*> gmhist_2D;
  vector<TH2*> gmhist_bUp_2D;
  vector<TH2*> gmhist_bDw_2D;
  vector<TH2*> gmhist_metJetUp_2D;
  vector<TH2*> gmhist_metJetDw_2D;
  vector<TH2*> gmhist_metResUp_2D;
  vector<TH2*> gmhist_metResDw_2D;
  vector<TH2*> gmhist_metUncUp_2D;
  vector<TH2*> gmhist_metUncDw_2D; 
  vector<TH2*> vlhist_2D; 
  vector<TH2*> vlhist_bUp_2D;
  vector<TH2*> vlhist_bDw_2D;
  vector<TH2*> vlhist_metJetUp_2D;
  vector<TH2*> vlhist_metJetDw_2D;
  vector<TH2*> vlhist_metResUp_2D;
  vector<TH2*> vlhist_metResDw_2D;
  vector<TH2*> vlhist_metUncUp_2D;
  vector<TH2*> vlhist_metUncDw_2D;
  vector<TH2*> vllhist_2D;
  vector<TH2*> vllhist_bUp_2D;
  vector<TH2*> vllhist_bDw_2D;
  vector<TH2*> vllhist_metJetUp_2D;
  vector<TH2*> vllhist_metJetDw_2D;
  vector<TH2*> vllhist_metResUp_2D;
  vector<TH2*> vllhist_metResDw_2D;
  vector<TH2*> vllhist_metUncUp_2D;
  vector<TH2*> vllhist_metUncDw_2D;

  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;


    TH2F* dthist_temp  = new TH2F((string("datahist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* tthist_temp  = new TH2F((string("tbkghist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* dbhist_temp  = new TH2F((string("dbkghist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* gmhist_temp  = new TH2F((string("gbkghist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* qchist_temp  = new TH2F((string("qbkghist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* vlhist_temp  = new TH2F((string("vlbkghist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* vllhist_temp = new TH2F((string("vllbkghist")+suffix+"_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    dbhist_2D.push_back(dynamic_cast<TH2*>(dbhist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    qchist_2D.push_back(dynamic_cast<TH2*>(qchist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vllhist_2D.push_back(dynamic_cast<TH2*>(vllhist_temp));

    if(doShapeSystematics){

      TH2F* tthist_bUp_temp = new TH2F((string("tbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_bDw_temp = new TH2F((string("tbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_metJetUp_temp = new TH2F((string("tbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_metJetDw_temp = new TH2F((string("tbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_metResUp_temp = new TH2F((string("tbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_metResDw_temp = new TH2F((string("tbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_metUncUp_temp = new TH2F((string("tbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* tthist_metUncDw_temp = new TH2F((string("tbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

      tthist_bUp_2D.push_back(dynamic_cast<TH2*>(tthist_bUp_temp));
      tthist_bDw_2D.push_back(dynamic_cast<TH2*>(tthist_bDw_temp));
      tthist_metJetUp_2D.push_back(dynamic_cast<TH2*>(tthist_metJetUp_temp));
      tthist_metJetDw_2D.push_back(dynamic_cast<TH2*>(tthist_metJetDw_temp));
      tthist_metResUp_2D.push_back(dynamic_cast<TH2*>(tthist_metResUp_temp));
      tthist_metResDw_2D.push_back(dynamic_cast<TH2*>(tthist_metResDw_temp));
      tthist_metUncUp_2D.push_back(dynamic_cast<TH2*>(tthist_metUncUp_temp));
      tthist_metUncDw_2D.push_back(dynamic_cast<TH2*>(tthist_metUncDw_temp));

      TH2F* dbhist_bUp_temp = new TH2F((string("dbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_bDw_temp = new TH2F((string("dbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metJetUp_temp = new TH2F((string("dbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metJetDw_temp = new TH2F((string("dbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metResUp_temp = new TH2F((string("dbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metResDw_temp = new TH2F((string("dbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metUncUp_temp = new TH2F((string("dbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metUncDw_temp = new TH2F((string("dbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

      dbhist_bUp_2D.push_back(dynamic_cast<TH2*>(dbhist_bUp_temp));
      dbhist_bDw_2D.push_back(dynamic_cast<TH2*>(dbhist_bDw_temp));
      dbhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetUp_temp));
      dbhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetDw_temp));
      dbhist_metResUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metResUp_temp));
      dbhist_metResDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metResDw_temp));
      dbhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncUp_temp));
      dbhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncDw_temp));

      TH2F* gmhist_bUp_temp = new TH2F((string("gbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_bDw_temp = new TH2F((string("gbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metJetUp_temp = new TH2F((string("gbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metJetDw_temp = new TH2F((string("gbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metResUp_temp = new TH2F((string("gbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metResDw_temp = new TH2F((string("gbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metUncUp_temp = new TH2F((string("gbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metUncDw_temp = new TH2F((string("gbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

      gmhist_bUp_2D.push_back(dynamic_cast<TH2*>(gmhist_bUp_temp));
      gmhist_bDw_2D.push_back(dynamic_cast<TH2*>(gmhist_bDw_temp));
      gmhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetUp_temp));
      gmhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetDw_temp));
      gmhist_metResUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metResUp_temp));
      gmhist_metResDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metResDw_temp));
      gmhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncUp_temp));
      gmhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncDw_temp));

    }

    if((sample == 1 or sample == 3) and doShapeSystematics){
      TH2F* vlhist_bUp_temp = new TH2F((string("vlbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_bDw_temp = new TH2F((string("vlbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metJetUp_temp = new TH2F((string("vlbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metJetDw_temp = new TH2F((string("vlbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metResUp_temp = new TH2F((string("vlbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metResDw_temp = new TH2F((string("vlbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metUncUp_temp = new TH2F((string("vlbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metUncDw_temp = new TH2F((string("vlbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

      vlhist_bUp_2D.push_back(dynamic_cast<TH2*>(vlhist_bUp_temp));
      vlhist_bDw_2D.push_back(dynamic_cast<TH2*>(vlhist_bDw_temp));
      vlhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetUp_temp));
      vlhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetDw_temp));
      vlhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metResUp_temp));
      vlhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metResDw_temp));
      vlhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncUp_temp));
      vlhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncDw_temp));
      
    }
    else if((sample == 2 or sample == 4) and doShapeSystematics){

      TH2F* vllhist_bUp_temp = new TH2F((string("vllbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_bDw_temp = new TH2F((string("vllbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metJetUp_temp = new TH2F((string("vllbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metJetDw_temp = new TH2F((string("vllbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metResUp_temp = new TH2F((string("vllbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metResDw_temp = new TH2F((string("vllbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metUncUp_temp = new TH2F((string("vllbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metUncDw_temp = new TH2F((string("vllbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

      vllhist_bUp_2D.push_back(dynamic_cast<TH2*>(vllhist_bUp_temp));
      vllhist_bDw_2D.push_back(dynamic_cast<TH2*>(vllhist_bDw_temp));
      vllhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetUp_temp));
      vllhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetDw_temp));
      vllhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metResUp_temp));
      vllhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metResDw_temp));
      vllhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncUp_temp));
      vllhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncDw_temp));
     
    }
  }

  TChain* dttree = new TChain("tree");
  dttree->Add(dtfile->GetName());
  if(dtfile_2)
    dttree->Add(dtfile_2->GetName());
  
  TTree* vltree  = (TTree*)vlfile->Get("tree/tree");
  TTree* vlltree = (TTree*)vllfile->Get("tree/tree");
  TTree* tttree  = (TTree*)ttfile->Get("tree/tree");
  TTree* dbtree  = (TTree*)dbfile->Get("tree/tree");
  TTree* gmtree  = (TTree*)gmfile->Get("tree/tree");
  TTree* qctree  = (TTree*)qcfile->Get("tree/tree");


  TFile kffile(kfactorFile.c_str());  

  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  if(not znlohist or not zlohist or not zewkhist){
    znlohist = (TH1*) kffile.Get("znlo012/znlo012_nominal");
    zlohist  = (TH1*) kffile.Get("zlo/zlo_nominal");
    zewkhist = (TH1*) kffile.Get("z_ewkcorr/z_ewkcorr");
    znlohist->Divide(zlohist);
  }

  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  if(not wnlohist or not wlohist or not wewkhist){
    wnlohist = (TH1*) kffile.Get("wnlo012/wnlo012_nominal");
    wlohist  = (TH1*) kffile.Get("wlo/wlo_nominal");
    wewkhist = (TH1*) kffile.Get("w_ewkcorr/w_ewkcorr");
    wnlohist->Divide(wlohist);
  }

  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  if(not anlohist or not alohist or not aewkhist){
    anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
    alohist  = (TH1*)kffile.Get("alo/alo_nominal");
    aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");
    anlohist->Divide(alohist);
  }

  vector<TH1*> ehists;
  vector<TH1*> vlhists;
  vector<TH1*> ahists;
  vector<TH1*> vllhists;
  // apply NLO QCD and EWK corrections for Zll and Wlnu                                                                                                                      
  ahists.push_back(anlohist);
  ahists.push_back(aewkhist);
  vlhists.push_back(wnlohist);
  vlhists.push_back(wewkhist);
  vllhists.push_back(znlohist);
  vllhists.push_back(zewkhist);

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  int indexQGL_W = 0;
  int indexQGL_Z = 0;
  int indexQGL_G = 0;
  int indexQGL_T = 0;

  if(applyQGLReweight){
    indexQGL_G = 3;
    indexQGL_W = 2;
    indexQGL_Z = 1;
    indexQGL_T = 4;
  }

  cout<<"lepton+jets control region --> W+jets"<<endl;
  makehist4(vltree, vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "", false, true,0,isHiggsInvisible);
  cout<<"lepton+jets control region --> Z+jets"<<endl;
  makehist4(vlltree,vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists, "", false, true,0,isHiggsInvisible);
  cout<<"lepton+jets control region --> top"<<endl;
  makehist4(tttree, tthist,  tthist_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,   "", true, true,0,isHiggsInvisible);
  cout<<"lepton+jets control region --> Diboson"<<endl;
  makehist4(dbtree, dbhist,  dbhist_2D,  true,  sample, category, isWJet, 1.00, lumi, 0, ehists, "", false, true,0,isHiggsInvisible);
  cout<<"lepton+jets control region --> gamma+jets"<<endl;
  makehist4(gmtree, gmhist,  gmhist_2D,  true,  sample, category, false, 1.00, lumi, indexQGL_G, ahists, "", false, true,0,isHiggsInvisible);
  cout<<"lepton+jets control region --> QCD"<<endl;
  makehist4(qctree, qchist,  qchist_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists, "", false,true,0,isHiggsInvisible);

  if(doShapeSystematics and (sample == 1 or sample == 3)){
    cout<<"lepton +jets region --> systematics for W+jets"<<endl;
    makehist4(vltree, vlhist_bUp,  vlhist_bUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "btagUp", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_bDw,  vlhist_bDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "btagDown", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_metJetUp,  vlhist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "jesUp", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_metJetDw,  vlhist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "jesDw", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_metResUp,  vlhist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "jerUp", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_metResDw,  vlhist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "jerDw", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_metUncUp,  vlhist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "uncUp", false, true,0,isHiggsInvisible);
    makehist4(vltree, vlhist_metUncDw,  vlhist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_W, vlhists,  "uncDw", false, true,0,isHiggsInvisible);
  }
  else if(doShapeSystematics and (sample == 2 or sample == 4)){
    cout<<"lepton +jets region --> systematics for Z+jets"<<endl;
    makehist4(vlltree, vllhist_bUp,  vllhist_bUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "btagUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_bDw,  vllhist_bDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "btagDown", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_metJetUp,  vllhist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "jesUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_metJetDw,  vllhist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "jesDw", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_metResUp,  vllhist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "jerUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_metResDw,  vllhist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "jerDw", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_metUncUp,  vllhist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "uncUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree, vllhist_metUncDw,  vllhist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_Z, vllhists,  "uncDw", false, true,0,isHiggsInvisible);
  }

  if(doShapeSystematics){

    cout<<"lepton +jets region --> systematics for top"<<endl;
    makehist4(tttree, tthist_bUp,  tthist_bUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "btagUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_bDw,  tthist_bDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "btagDown", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metJetUp,  tthist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "jesUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metJetDw,  tthist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "jesDw", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metResUp,  tthist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "jerUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metResDw,  tthist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "jerDw", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metUncUp,  tthist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "uncUp", true, true,0,isHiggsInvisible);
    makehist4(tttree, tthist_metUncDw,  tthist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, indexQGL_T, ehists,  "uncDw", true, true,0,isHiggsInvisible);

    cout<<"lepton +jets region --> systematics for di-boson"<<endl;
    makehist4(dbtree, dbhist_bUp,  dbhist_bUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "btagUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_bDw,  dbhist_bDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "btagDown", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_metJetUp,  dbhist_metJetUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "jesUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_metJetDw,  dbhist_metJetDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "jesDw", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_metResUp,  dbhist_metResUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "jerUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_metResDw,  dbhist_metResDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "jerDw", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_metUncUp,  dbhist_metUncUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "uncUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree, dbhist_metUncDw,  dbhist_metUncDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists,  "uncDw", false, true,0,isHiggsInvisible);

    cout<<"lepton +jets region --> systematics for gamma+jets"<<endl;
    makehist4(gmtree, gmhist_bUp,  gmhist_bUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "btagUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_bDw,  gmhist_bDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "btagDown", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metJetUp,  gmhist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "jesUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metJetDw,  gmhist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "jesDw", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metResUp,  gmhist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "jerUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metResDw,  gmhist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "jerDw", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metUncUp,  gmhist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "uncUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree, gmhist_metUncDw,  gmhist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists,  "uncDw", false, true,0,isHiggsInvisible);


  }
  
  cout<<"lepton+jets control region --> Data"<<endl;
  makehist4(dttree, dthist, dthist_2D,   false, sample, category, false,  1.00, lumi, 0, ehists, "", false, true,0,isHiggsInvisible);

  //smooth
  for(auto hist : tthist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dbhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : gmhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vlhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vllhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : qchist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  for(auto hist : tthist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dbhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : gmhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vlhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vllhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : qchist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}


  if(doShapeSystematics){

    for(auto hist : tthist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : dbhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vlhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vllhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : tthist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : tthist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : dbhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vlhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vllhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
   
  }
  
  //
  outfile->cd();
  for(auto hist :  dthist) hist->Write();
  for(auto hist :  tthist) hist->Write();
  for(auto hist :  dbhist) hist->Write();
  for(auto hist :  gmhist) hist->Write();
  for(auto hist :  qchist) hist->Write();
  for(auto hist :  vlhist) hist->Write();
  for(auto hist :  vllhist) hist->Write();

  if(doShapeSystematics){
    
    for(auto hist :  dbhist_bUp) hist->Write();
    for(auto hist :  dbhist_bDw) hist->Write();
    for(auto hist :  dbhist_metJetUp) hist->Write();
    for(auto hist :  dbhist_metJetDw) hist->Write();
    for(auto hist :  dbhist_metResUp) hist->Write();
    for(auto hist :  dbhist_metResDw) hist->Write();
    for(auto hist :  dbhist_metUncUp) hist->Write();
    for(auto hist :  dbhist_metUncDw) hist->Write();

    for(auto hist :  gmhist_bUp) hist->Write();
    for(auto hist :  gmhist_bDw) hist->Write();
    for(auto hist :  gmhist_metJetUp) hist->Write();
    for(auto hist :  gmhist_metJetDw) hist->Write();
    for(auto hist :  gmhist_metResUp) hist->Write();
    for(auto hist :  gmhist_metResDw) hist->Write();
    for(auto hist :  gmhist_metUncUp) hist->Write();
    for(auto hist :  gmhist_metUncDw) hist->Write();
 
    for(auto hist :  tthist_bUp) hist->Write();
    for(auto hist :  tthist_bDw) hist->Write();
    for(auto hist :  tthist_metJetUp) hist->Write();
    for(auto hist :  tthist_metJetDw) hist->Write();
    for(auto hist :  tthist_metResUp) hist->Write();
    for(auto hist :  tthist_metResDw) hist->Write();
    for(auto hist :  tthist_metUncUp) hist->Write();
    for(auto hist :  tthist_metUncDw) hist->Write();

    for(auto hist :  vlhist_bUp) hist->Write();
    for(auto hist :  vlhist_bDw) hist->Write();
    for(auto hist :  vlhist_metJetUp) hist->Write();
    for(auto hist :  vlhist_metJetDw) hist->Write();
    for(auto hist :  vlhist_metResUp) hist->Write();
    for(auto hist :  vlhist_metResDw) hist->Write();
    for(auto hist :  vlhist_metUncUp) hist->Write();
    for(auto hist :  vlhist_metUncDw) hist->Write();

    for(auto hist :  vllhist_bUp) hist->Write();
    for(auto hist :  vllhist_bDw) hist->Write();
    for(auto hist :  vllhist_metJetUp) hist->Write();
    for(auto hist :  vllhist_metJetDw) hist->Write();
    for(auto hist :  vllhist_metResUp) hist->Write();
    for(auto hist :  vllhist_metResDw) hist->Write();
    for(auto hist :  vllhist_metUncUp) hist->Write();
    for(auto hist :  vllhist_metUncDw) hist->Write();
  }

  for(auto hist_2D : vlhist_2D) hist_2D->Write();
  for(auto hist_2D : vllhist_2D) hist_2D->Write();
  for(auto hist_2D : tthist_2D) hist_2D->Write();
  for(auto hist_2D : dbhist_2D) hist_2D->Write();
  for(auto hist_2D : gmhist_2D) hist_2D->Write();
  for(auto hist_2D : qchist_2D) hist_2D->Write();
  for(auto hist_2D : dthist_2D) hist_2D->Write();

  for(auto hist_2D : vlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vllhist_2D){TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : tthist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dbhist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : gmhist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : qchist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dthist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

  if(doShapeSystematics){
    
    for(auto hist :  dbhist_bUp_2D) hist->Write();
    for(auto hist :  dbhist_bDw_2D) hist->Write();
    for(auto hist :  dbhist_metJetUp_2D) hist->Write();
    for(auto hist :  dbhist_metJetDw_2D) hist->Write();
    for(auto hist :  dbhist_metResUp_2D) hist->Write();
    for(auto hist :  dbhist_metResDw_2D) hist->Write();
    for(auto hist :  dbhist_metUncUp_2D) hist->Write();
    for(auto hist :  dbhist_metUncDw_2D) hist->Write();

    for(auto hist : dbhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

    for(auto hist :  gmhist_bUp_2D) hist->Write();
    for(auto hist :  gmhist_bDw_2D) hist->Write();
    for(auto hist :  gmhist_metJetUp_2D) hist->Write();
    for(auto hist :  gmhist_metJetDw_2D) hist->Write();
    for(auto hist :  gmhist_metResUp_2D) hist->Write();
    for(auto hist :  gmhist_metResDw_2D) hist->Write();
    for(auto hist :  gmhist_metUncUp_2D) hist->Write();
    for(auto hist :  gmhist_metUncDw_2D) hist->Write();

    for(auto hist : gmhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
 
    for(auto hist :  tthist_bUp_2D) hist->Write();
    for(auto hist :  tthist_bDw_2D) hist->Write();
    for(auto hist :  tthist_metJetUp_2D) hist->Write();
    for(auto hist :  tthist_metJetDw_2D) hist->Write();
    for(auto hist :  tthist_metResUp_2D) hist->Write();
    for(auto hist :  tthist_metResDw_2D) hist->Write();
    for(auto hist :  tthist_metUncUp_2D) hist->Write();
    for(auto hist :  tthist_metUncDw_2D) hist->Write();

    for(auto hist : tthist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

    for(auto hist :  vlhist_bUp_2D) hist->Write();
    for(auto hist :  vlhist_bDw_2D) hist->Write();
    for(auto hist :  vlhist_metJetUp_2D) hist->Write();
    for(auto hist :  vlhist_metJetDw_2D) hist->Write();
    for(auto hist :  vlhist_metResUp_2D) hist->Write();
    for(auto hist :  vlhist_metResDw_2D) hist->Write();
    for(auto hist :  vlhist_metUncUp_2D) hist->Write();
    for(auto hist :  vlhist_metUncDw_2D) hist->Write();

    for(auto hist : vlhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

    for(auto hist :  vllhist_bUp_2D) hist->Write();
    for(auto hist :  vllhist_bDw_2D) hist->Write();
    for(auto hist :  vllhist_metJetUp_2D) hist->Write();
    for(auto hist :  vllhist_metJetDw_2D) hist->Write();
    for(auto hist :  vllhist_metResUp_2D) hist->Write();
    for(auto hist :  vllhist_metResDw_2D) hist->Write();
    for(auto hist :  vllhist_metUncUp_2D) hist->Write();
    for(auto hist :  vllhist_metUncDw_2D) hist->Write();

    for(auto hist : vllhist_bUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_bDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

  }


  dtfile->Close();
  vlfile->Close();
  vllfile->Close();
  ttfile->Close();
  dbfile->Close();
  gmfile->Close();
  qcfile->Close();

  dthist.clear();
  tthist.clear();
  tthist_bUp.clear();
  tthist_bDw.clear();
  tthist_metJetUp.clear();
  tthist_metJetDw.clear();
  tthist_metResUp.clear();
  tthist_metResDw.clear();
  tthist_metUncUp.clear();
  tthist_metUncDw.clear();
  qchist.clear();
  dbhist.clear();
  dbhist_bUp.clear();
  dbhist_bDw.clear();
  dbhist_metJetUp.clear();
  dbhist_metJetDw.clear();
  dbhist_metResUp.clear();
  dbhist_metResDw.clear();
  dbhist_metUncUp.clear();
  dbhist_metUncDw.clear();
  gmhist.clear();
  gmhist_bUp.clear();
  gmhist_bDw.clear();
  gmhist_metJetUp.clear();
  gmhist_metJetDw.clear();
  gmhist_metResUp.clear();
  gmhist_metResDw.clear();
  gmhist_metUncUp.clear();
  gmhist_metUncDw.clear();
  vlhist.clear();
  vlhist_bUp.clear();
  vlhist_bDw.clear();
  vlhist_metJetUp.clear();
  vlhist_metJetDw.clear();
  vlhist_metResUp.clear();
  vlhist_metResDw.clear();
  vlhist_metUncUp.clear();
  vlhist_metUncDw.clear();
  vllhist.clear();
  vllhist_bUp.clear();
  vllhist_bDw.clear();
  vllhist_metJetUp.clear();
  vllhist_metJetDw.clear();
  vllhist_metResUp.clear();
  vllhist_metResDw.clear();
  vllhist_metUncUp.clear();
  vllhist_metUncDw.clear();


  dthist_2D.clear();
  tthist_2D.clear();
  tthist_bUp_2D.clear();
  tthist_bDw_2D.clear();
  tthist_metJetUp_2D.clear();
  tthist_metJetDw_2D.clear();
  tthist_metResUp_2D.clear();
  tthist_metResDw_2D.clear();
  tthist_metUncUp_2D.clear();
  tthist_metUncDw_2D.clear();
  qchist_2D.clear();
  dbhist_2D.clear();
  dbhist_bUp_2D.clear();
  dbhist_bDw_2D.clear();
  dbhist_metJetUp_2D.clear();
  dbhist_metJetDw_2D.clear();
  dbhist_metResUp_2D.clear();
  dbhist_metResDw_2D.clear();
  dbhist_metUncUp_2D.clear();
  dbhist_metUncDw_2D.clear();
  gmhist_2D.clear();
  gmhist_bUp_2D.clear();
  gmhist_bDw_2D.clear();
  gmhist_metJetUp_2D.clear();
  gmhist_metJetDw_2D.clear();
  gmhist_metResUp_2D.clear();
  gmhist_metResDw_2D.clear();
  gmhist_metUncUp_2D.clear();
  gmhist_metUncDw_2D.clear();
  vlhist_2D.clear();
  vlhist_bUp_2D.clear();
  vlhist_bDw_2D.clear();
  vlhist_metJetUp_2D.clear();
  vlhist_metJetDw_2D.clear();
  vlhist_metResUp_2D.clear();
  vlhist_metResDw_2D.clear();
  vlhist_metUncUp_2D.clear();
  vlhist_metUncDw_2D.clear();
  vllhist_2D.clear();
  vllhist_bUp_2D.clear();
  vllhist_bDw_2D.clear();
  vllhist_metJetUp_2D.clear();
  vllhist_metJetDw_2D.clear();
  vllhist_metResUp_2D.clear();
  vllhist_metResDw_2D.clear();
  vllhist_metUncUp_2D.clear();
  vllhist_metUncDw_2D.clear();

  cout << "Templates for the lepton control region computed ..." << endl;
}


//build template for tt
void topdatamchist(TFile* outfile, int sample, int category, vector<string> observables, vector<string> observables_2D,
		   double lumi = 2.24, bool applyQGLReweight = false,
		   bool makeResonantSelection = false, bool doShapeSystematics = false, bool isHiggsInvisible = false) {

  if (sample != 7 && sample != 8) return;

  TFile* ttfile      = NULL;
  TFile* ttfile_alt  = NULL;
  TFile* dbfile  = NULL;
  TFile* gmfile  = NULL;
  TFile* qcfile  = NULL;
  TFile* vlfile  = NULL;
  TFile* vllfile = NULL;
  TFile* dtfile  = NULL;

  string suffix;

  if(sample == 7)
    suffix = "topmu";
  else if(sample == 8)
    suffix = "topel";

  vllfile = TFile::Open((baseInputTreePath+"DYJets/topfilter/top_tree_DYJetsToLL_M-50.root").c_str());
  vlfile = TFile::Open((baseInputTreePath+"WJets/topfilter/top_tree_WJetsToLNu.root").c_str());
  qcfile = TFile::Open((baseInputTreePath+"QCD/topfilter/top_tree_QCD.root").c_str());
  dbfile = TFile::Open((baseInputTreePath+"DiBoson/topfilter/top_tree_DiBoson.root").c_str());
  gmfile = TFile::Open((baseInputTreePath+"DiBoson/topfilter/top_tree_DiBoson.root").c_str());
  ttfile = TFile::Open((baseInputTreePath+"Top/topfilter/top_tree_Top_amc.root").c_str());
  ttfile_alt = TFile::Open((baseInputTreePath+"Top/topfilter/top_tree_Top.root").c_str());

  if(sample == 7)
    dtfile = TFile::Open((baseInputTreePath+"MET/topfilter/top_tree_crab_MET-Run2015.root").c_str());
  else if(sample == 8)
    dtfile = TFile::Open((baseInputTreePath+"SingleElectron/topfilter/top_tree_crab_SingleEle-Run2015.root").c_str());

  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_alt;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> dbhist_bUp;
  vector<TH1*> dbhist_bDw;
  vector<TH1*> dbhist_metJetUp;
  vector<TH1*> dbhist_metJetDw;
  vector<TH1*> dbhist_metResUp;
  vector<TH1*> dbhist_metResDw;
  vector<TH1*> dbhist_metUncUp;
  vector<TH1*> dbhist_metUncDw;

  vector<TH1*> gmhist;
  vector<TH1*> gmhist_bUp;
  vector<TH1*> gmhist_bDw;
  vector<TH1*> gmhist_metJetUp;
  vector<TH1*> gmhist_metJetDw;
  vector<TH1*> gmhist_metResUp;
  vector<TH1*> gmhist_metResDw;
  vector<TH1*> gmhist_metUncUp;
  vector<TH1*> gmhist_metUncDw;

  vector<TH1*> vlhist;
  vector<TH1*> vlhist_bUp;
  vector<TH1*> vlhist_bDw;
  vector<TH1*> vlhist_metJetUp;
  vector<TH1*> vlhist_metJetDw;
  vector<TH1*> vlhist_metResUp;
  vector<TH1*> vlhist_metResDw;
  vector<TH1*> vlhist_metUncUp;
  vector<TH1*> vlhist_metUncDw;
  vector<TH1*> vllhist;
  vector<TH1*> vllhist_bUp;
  vector<TH1*> vllhist_bDw;
  vector<TH1*> vllhist_metJetUp;
  vector<TH1*> vllhist_metJetDw;
  vector<TH1*> vllhist_metResUp;
  vector<TH1*> vllhist_metResDw;
  vector<TH1*> vllhist_metUncUp;
  vector<TH1*> vllhist_metUncDw;

  vector<TH1*> tthist_matched;
  vector<TH1*> tthist_matched_alt;
  vector<TH1*> tthist_unmatched;
  vector<TH1*> tthist_unmatched_alt;

  vector<float> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "",int(bins.size()-1), &bins[0]);
    TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "",int(bins.size()-1), &bins[0]);
    TH1F* gmhist_temp = new TH1F((string("gbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
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

    if(doShapeSystematics){
      
      TH1F* dbhist_bUp_temp = new TH1F((string("dbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_bDw_temp = new TH1F((string("dbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metJetUp_temp = new TH1F((string("dbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metJetDw_temp = new TH1F((string("dbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metResUp_temp = new TH1F((string("dbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metResDw_temp = new TH1F((string("dbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metUncUp_temp = new TH1F((string("dbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_metUncDw_temp = new TH1F((string("dbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      
      dbhist_bUp.push_back(dynamic_cast<TH1*>(dbhist_bUp_temp));
      dbhist_bDw.push_back(dynamic_cast<TH1*>(dbhist_bDw_temp));
      dbhist_metJetUp.push_back(dynamic_cast<TH1*>(dbhist_metJetUp_temp));
      dbhist_metJetDw.push_back(dynamic_cast<TH1*>(dbhist_metJetDw_temp));
      dbhist_metResUp.push_back(dynamic_cast<TH1*>(dbhist_metResUp_temp));
      dbhist_metResDw.push_back(dynamic_cast<TH1*>(dbhist_metResDw_temp));
      dbhist_metUncUp.push_back(dynamic_cast<TH1*>(dbhist_metUncUp_temp));
      dbhist_metUncDw.push_back(dynamic_cast<TH1*>(dbhist_metUncDw_temp));


      TH1F* gmhist_bUp_temp = new TH1F((string("gbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_bDw_temp = new TH1F((string("gbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metJetUp_temp = new TH1F((string("gbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metJetDw_temp = new TH1F((string("gbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metResUp_temp = new TH1F((string("gbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metResDw_temp = new TH1F((string("gbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metUncUp_temp = new TH1F((string("gbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* gmhist_metUncDw_temp = new TH1F((string("gbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      
      gmhist_bUp.push_back(dynamic_cast<TH1*>(gmhist_bUp_temp));
      gmhist_bDw.push_back(dynamic_cast<TH1*>(gmhist_bDw_temp));
      gmhist_metJetUp.push_back(dynamic_cast<TH1*>(gmhist_metJetUp_temp));
      gmhist_metJetDw.push_back(dynamic_cast<TH1*>(gmhist_metJetDw_temp));
      gmhist_metResUp.push_back(dynamic_cast<TH1*>(gmhist_metResUp_temp));
      gmhist_metResDw.push_back(dynamic_cast<TH1*>(gmhist_metResDw_temp));
      gmhist_metUncUp.push_back(dynamic_cast<TH1*>(gmhist_metUncUp_temp));
      gmhist_metUncDw.push_back(dynamic_cast<TH1*>(gmhist_metUncDw_temp));

      
      TH1F* vlhist_bUp_temp = new TH1F((string("vlbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_bDw_temp = new TH1F((string("vlbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metJetUp_temp = new TH1F((string("vlbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metJetDw_temp = new TH1F((string("vlbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metResUp_temp = new TH1F((string("vlbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metResDw_temp = new TH1F((string("vlbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metUncUp_temp = new TH1F((string("vlbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_metUncDw_temp = new TH1F((string("vlbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      vlhist_bUp.push_back(dynamic_cast<TH1*>(vlhist_bUp_temp));
      vlhist_bDw.push_back(dynamic_cast<TH1*>(vlhist_bDw_temp));
      vlhist_metJetUp.push_back(dynamic_cast<TH1*>(vlhist_metJetUp_temp));
      vlhist_metJetDw.push_back(dynamic_cast<TH1*>(vlhist_metJetDw_temp));
      vlhist_metResUp.push_back(dynamic_cast<TH1*>(vlhist_metResUp_temp));
      vlhist_metResDw.push_back(dynamic_cast<TH1*>(vlhist_metResDw_temp));
      vlhist_metUncUp.push_back(dynamic_cast<TH1*>(vlhist_metUncUp_temp));
      vlhist_metUncDw.push_back(dynamic_cast<TH1*>(vlhist_metUncDw_temp));

      TH1F* vllhist_bUp_temp = new TH1F((string("vllbkghist")+suffix+"_bUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_bDw_temp = new TH1F((string("vllbkghist")+suffix+"_bDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metJetUp_temp = new TH1F((string("vllbkghist")+suffix+"_metJetUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metJetDw_temp = new TH1F((string("vllbkghist")+suffix+"_metJetDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metResUp_temp = new TH1F((string("vllbkghist")+suffix+"_metResUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metResDw_temp = new TH1F((string("vllbkghist")+suffix+"_metResDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metUncUp_temp = new TH1F((string("vllbkghist")+suffix+"_metUncUp_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_metUncDw_temp = new TH1F((string("vllbkghist")+suffix+"_metUncDw_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
 
      vllhist_bUp.push_back(dynamic_cast<TH1*>(vllhist_bUp_temp));
      vllhist_bDw.push_back(dynamic_cast<TH1*>(vllhist_bDw_temp));
      vllhist_metJetUp.push_back(dynamic_cast<TH1*>(vllhist_metJetUp_temp));
      vllhist_metJetDw.push_back(dynamic_cast<TH1*>(vllhist_metJetDw_temp));
      vllhist_metResUp.push_back(dynamic_cast<TH1*>(vllhist_metResUp_temp));
      vllhist_metResDw.push_back(dynamic_cast<TH1*>(vllhist_metResDw_temp));
      vllhist_metUncUp.push_back(dynamic_cast<TH1*>(vllhist_metUncUp_temp));
      vllhist_metUncDw.push_back(dynamic_cast<TH1*>(vllhist_metUncDw_temp));
      
    }

  }

  vector<TH2*> dthist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_alt_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> dbhist_bUp_2D;
  vector<TH2*> dbhist_bDw_2D;
  vector<TH2*> dbhist_metJetUp_2D;
  vector<TH2*> dbhist_metJetDw_2D;
  vector<TH2*> dbhist_metResUp_2D;
  vector<TH2*> dbhist_metResDw_2D;
  vector<TH2*> dbhist_metUncUp_2D;
  vector<TH2*> dbhist_metUncDw_2D;
  vector<TH2*> gmhist_2D;
  vector<TH2*> gmhist_bUp_2D;
  vector<TH2*> gmhist_bDw_2D;
  vector<TH2*> gmhist_metJetUp_2D;
  vector<TH2*> gmhist_metJetDw_2D;
  vector<TH2*> gmhist_metResUp_2D;
  vector<TH2*> gmhist_metResDw_2D;
  vector<TH2*> gmhist_metUncUp_2D;
  vector<TH2*> gmhist_metUncDw_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vlhist_bUp_2D;
  vector<TH2*> vlhist_bDw_2D;
  vector<TH2*> vlhist_metJetUp_2D;
  vector<TH2*> vlhist_metJetDw_2D;
  vector<TH2*> vlhist_metResUp_2D;
  vector<TH2*> vlhist_metResDw_2D;
  vector<TH2*> vlhist_metUncUp_2D;
  vector<TH2*> vlhist_metUncDw_2D;
  vector<TH2*> vllhist_2D;
  vector<TH2*> vllhist_bUp_2D;
  vector<TH2*> vllhist_bDw_2D;
  vector<TH2*> vllhist_metJetUp_2D;
  vector<TH2*> vllhist_metJetDw_2D;
  vector<TH2*> vllhist_metResUp_2D;
  vector<TH2*> vllhist_metResDw_2D;
  vector<TH2*> vllhist_metUncUp_2D;
  vector<TH2*> vllhist_metUncDw_2D;

  vector<TH2*> tthist_matched_2D;
  vector<TH2*> tthist_matched_alt_2D;
  vector<TH2*> tthist_unmatched_2D;
  vector<TH2*> tthist_unmatched_alt_2D;

  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable --> please define it"<<endl;

    TH2F* dthist_temp     = new TH2F((string("datahist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* dbhist_temp     = new TH2F((string("dbkghist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* gmhist_temp     = new TH2F((string("gbkghist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* qchist_temp     = new TH2F((string("qbkghist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* vlhist_temp     = new TH2F((string("vlbkghist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* vllhist_temp    = new TH2F((string("vllbkghist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* tthist_temp     = new TH2F((string("tbkghist")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
    TH2F* tthist_alt_temp = new TH2F((string("tbkghist_alt")+suffix+"_"+obs+"_2D").c_str(), "",int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
  
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    tthist_alt_2D.push_back(dynamic_cast<TH2*>(tthist_alt_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));
    dbhist_2D.push_back(dynamic_cast<TH2*>(dbhist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    qchist_2D.push_back(dynamic_cast<TH2*>(qchist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vllhist_2D.push_back(dynamic_cast<TH2*>(vllhist_temp));

    if(doShapeSystematics){
      
      TH2F* dbhist_bUp_temp = new TH2F((string("dbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_bDw_temp = new TH2F((string("dbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metJetUp_temp = new TH2F((string("dbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metJetDw_temp = new TH2F((string("dbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metResUp_temp = new TH2F((string("dbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metResDw_temp = new TH2F((string("dbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metUncUp_temp = new TH2F((string("dbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* dbhist_metUncDw_temp = new TH2F((string("dbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      
      dbhist_bUp_2D.push_back(dynamic_cast<TH2*>(dbhist_bUp_temp));
      dbhist_bDw_2D.push_back(dynamic_cast<TH2*>(dbhist_bDw_temp));
      dbhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetUp_temp));
      dbhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetDw_temp));
      dbhist_metResUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metResUp_temp));
      dbhist_metResDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metResDw_temp));
      dbhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncUp_temp));
      dbhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncDw_temp));


      TH2F* gmhist_bUp_temp = new TH2F((string("gbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_bDw_temp = new TH2F((string("gbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metJetUp_temp = new TH2F((string("gbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metJetDw_temp = new TH2F((string("gbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metResUp_temp = new TH2F((string("gbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metResDw_temp = new TH2F((string("gbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metUncUp_temp = new TH2F((string("gbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* gmhist_metUncDw_temp = new TH2F((string("gbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      
      gmhist_bUp_2D.push_back(dynamic_cast<TH2*>(gmhist_bUp_temp));
      gmhist_bDw_2D.push_back(dynamic_cast<TH2*>(gmhist_bDw_temp));
      gmhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetUp_temp));
      gmhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetDw_temp));
      gmhist_metResUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metResUp_temp));
      gmhist_metResDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metResDw_temp));
      gmhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncUp_temp));
      gmhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncDw_temp));

      
      TH2F* vlhist_bUp_temp = new TH2F((string("vlbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_bDw_temp = new TH2F((string("vlbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metJetUp_temp = new TH2F((string("vlbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metJetDw_temp = new TH2F((string("vlbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metResUp_temp = new TH2F((string("vlbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metResDw_temp = new TH2F((string("vlbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metUncUp_temp = new TH2F((string("vlbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vlhist_metUncDw_temp = new TH2F((string("vlbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);

      vlhist_bUp_2D.push_back(dynamic_cast<TH2*>(vlhist_bUp_temp));
      vlhist_bDw_2D.push_back(dynamic_cast<TH2*>(vlhist_bDw_temp));
      vlhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetUp_temp));
      vlhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetDw_temp));
      vlhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metResUp_temp));
      vlhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metResDw_temp));
      vlhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncUp_temp));
      vlhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncDw_temp));

      TH2F* vllhist_bUp_temp = new TH2F((string("vllbkghist")+suffix+"_bUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_bDw_temp = new TH2F((string("vllbkghist")+suffix+"_bDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metJetUp_temp = new TH2F((string("vllbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metJetDw_temp = new TH2F((string("vllbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metResUp_temp = new TH2F((string("vllbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metResDw_temp = new TH2F((string("vllbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metUncUp_temp = new TH2F((string("vllbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
      TH2F* vllhist_metUncDw_temp = new TH2F((string("vllbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(), "", int(bins.binX.size()-1), &bins.binX[0],int(bins.binY.size()-1), &bins.binY[0]);
 
      vllhist_bUp_2D.push_back(dynamic_cast<TH2*>(vllhist_bUp_temp));
      vllhist_bDw_2D.push_back(dynamic_cast<TH2*>(vllhist_bDw_temp));
      vllhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetUp_temp));
      vllhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetDw_temp));
      vllhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metResUp_temp));
      vllhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metResDw_temp));
      vllhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncUp_temp));
      vllhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncDw_temp));
      
    }
  }

  TTree* dttree  = (TTree*)dtfile->Get("tree");
  TTree* vltree  = (TTree*)vlfile->Get("tree/tree");
  TTree* vlltree = (TTree*)vllfile->Get("tree/tree");
  TTree* tttree  = (TTree*)ttfile->Get("tree/tree");
  TTree* tttree_alt = NULL;

  if(ttfile_alt)
    tttree_alt = (TTree*) ttfile_alt->Get("tree/tree");

  TTree* dbtree = (TTree*)dbfile->Get("tree/tree");
  TTree* gmtree = (TTree*)gmfile->Get("tree/tree");
  TTree* qctree = (TTree*)qcfile->Get("tree/tree");

  // k-factors
  TFile kffile(kfactorFile.c_str());  
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1* zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  if(not znlohist or not zlohist or not zewkhist){
    znlohist = (TH1*) kffile.Get("znlo012/znlo012_nominal");
    zlohist  = (TH1*) kffile.Get("zlo/zlo_nominal");
    zewkhist = (TH1*) kffile.Get("z_ewkcorr/z_ewkcorr");
    znlohist->Divide(zlohist);
  }

  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  if(not wnlohist or not wlohist or not wewkhist){
    wnlohist = (TH1*) kffile.Get("wnlo012/wnlo012_nominal");
    wlohist  = (TH1*) kffile.Get("wlo/wlo_nominal");
    wewkhist = (TH1*) kffile.Get("w_ewkcorr/w_ewkcorr");
    wnlohist->Divide(wlohist);
  }

  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  if(not anlohist or not alohist or not aewkhist){
    anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
    alohist  = (TH1*)kffile.Get("alo/alo_nominal");
    aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");
    anlohist->Divide(alohist);
  }

  vector<TH1*> ehists;
  vector<TH1*> vlhists;
  vector<TH1*> ahists;
  vector<TH1*> vllhists;
  // apply NLO QCD and EWK corrections for Zll and Wlnu                                                                                                                      
  ahists.push_back(anlohist);
  ahists.push_back(aewkhist);
  vllhists.push_back(znlohist);
  vllhists.push_back(zewkhist);
  vlhists.push_back(wnlohist);
  vlhists.push_back(wewkhist);

  // 
  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;


  int index_QGL_Z = 0;
  int index_QGL_G = 0;
  int index_QGL_W = 0;
  int index_QGL_T = 0;

  if(applyQGLReweight){
    index_QGL_Z = 1;
    index_QGL_W = 2;
    index_QGL_T = 4;
  }

  if(not makeResonantSelection){
    cout<<"top+jets control region --> Top"<<endl;
    makehist4(tttree,tthist,tthist_2D,true,  sample, category, false,  1.00, lumi, index_QGL_T, ehists, "", true, true,0,isHiggsInvisible);
    cout<<"top+jets control region --> Top alternative"<<endl;
    makehist4(tttree_alt,tthist_alt, tthist_alt_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_T, ehists, "", true, true,0,isHiggsInvisible);
  }
  else{
    cout<<"top+jets control region --> Top"<<endl;
    makehist4(tttree,tthist_matched,tthist_matched_2D,true, sample, category, false,  1.00, lumi, index_QGL_T, ehists, "", true, true,0,isHiggsInvisible,1);
    makehist4(tttree,tthist_unmatched,tthist_unmatched_2D,true,  sample, category, false,  1.00, lumi, index_QGL_T, ehists, "", true, true,0,isHiggsInvisible,2);
    cout<<"top+jets control region --> Top alternative"<<endl;
    makehist4(tttree_alt,tthist_matched_alt,tthist_matched_alt_2D,true,sample, category, false,  1.00, lumi, index_QGL_T, ehists, "", true, true,0,isHiggsInvisible,1);
    makehist4(tttree_alt,tthist_unmatched_alt,tthist_unmatched_alt_2D,true,sample,category, false,  1.00, lumi, index_QGL_T, ehists, "", true, true,0,isHiggsInvisible,2);
  }

  cout<<"top+jets control region --> W+jets"<<endl;
  makehist4(vltree,  vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "", false, true,0,isHiggsInvisible);
  cout<<"top+jets control region --> Z+jets"<<endl;
  makehist4(vlltree, vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "", false, true,0,isHiggsInvisible);
  cout<<"top+jets control region: Diboson"<<endl;
  makehist4(dbtree,  dbhist,  dbhist_2D,  true,  sample, category, isWJet, 1.00, lumi, 0, ehists, "", false, true,0,isHiggsInvisible);
  cout<<"top+jets control region: gamma+jets"<<endl;
  makehist4(gmtree,  gmhist,  gmhist_2D,  true,  sample, category, false, 1.00, lumi, index_QGL_G, ahists, "", false, true,0,isHiggsInvisible);
  cout<<"top+jets control region: QCD"<<endl;
  makehist4(qctree,  qchist,  qchist_2D,  true,  sample, category, false,  1.00, lumi, 0, ehists, "", false, true,0,isHiggsInvisible);

  if(doShapeSystematics){
    cout<<"top+jets control region -->  shape systematics for W+Jets"<<endl;
    makehist4(vltree,  vlhist_bUp,  vlhist_bUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "btagUp", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_bDw,  vlhist_bDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "btagDown", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_metJetUp,  vlhist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "jesUp", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_metJetDw,  vlhist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "jesDw", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_metResUp,  vlhist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "jerUp", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_metResDw,  vlhist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "jerDw", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_metUncUp,  vlhist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "uncUp", false, true,0,isHiggsInvisible);
    makehist4(vltree,  vlhist_metUncDw,  vlhist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_W, vlhists, "uncDw", false, true,0,isHiggsInvisible);

    cout<<"top+jets control region -->  shape systematics for Z+Jets"<<endl;
    makehist4(vlltree,  vllhist_bUp,  vllhist_bUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "btagUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_bDw,  vllhist_bDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "btagDown", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_metJetUp,  vllhist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "jesUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_metJetDw,  vllhist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "jesDw", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_metResUp,  vllhist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "jerUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_metResDw,  vllhist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "jerDw", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_metUncUp,  vllhist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "uncUp", false, true,0,isHiggsInvisible);
    makehist4(vlltree,  vllhist_metUncDw,  vllhist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, index_QGL_Z, vllhists, "uncDw", false, true,0,isHiggsInvisible);

    cout<<"top+jets control region -->  shape systematics for Dibosons"<<endl;
    makehist4(dbtree,  dbhist_bUp,  dbhist_bUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "btagUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_bDw,  dbhist_bDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "btagDown", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_metJetUp,  dbhist_metJetUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "jesUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_metJetDw,  dbhist_metJetDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "jesDw", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_metResUp,  dbhist_metResUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "jerUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_metResDw,  dbhist_metResDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "jerDw", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_metUncUp,  dbhist_metUncUp_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "uncUp", false, true,0,isHiggsInvisible);
    makehist4(dbtree,  dbhist_metUncDw,  dbhist_metUncDw_2D,  true,  sample, category, isWJet,  1.00, lumi, 0, ehists, "uncDw", false, true,0,isHiggsInvisible);

    cout<<"top+jets control region -->  shape systematics for gamma+jets"<<endl;
    makehist4(gmtree,  gmhist_bUp,  gmhist_bUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "btagUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_bDw,  gmhist_bDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "btagDown", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_metJetUp,  gmhist_metJetUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "jesUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_metJetDw,  gmhist_metJetDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "jesDw", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_metResUp,  gmhist_metResUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "jerUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_metResDw,  gmhist_metResDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "jerDw", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_metUncUp,  gmhist_metUncUp_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "uncUp", false, true,0,isHiggsInvisible);
    makehist4(gmtree,  gmhist_metUncDw,  gmhist_metUncDw_2D,  true,  sample, category, false,  1.00, lumi, 0, ahists, "uncDw", false, true,0,isHiggsInvisible);
  }


  cout<<"top+jets control region: Data"<<endl;
  makehist4(dttree,  dthist,  dthist_2D,  false, sample, category, false,  1.00, lumi, 0, ehists, "", false, true,0,isHiggsInvisible);
  
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


  //smooth
  for(auto hist : tthist){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dbhist){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : gmhist){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vlhist){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vllhist){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : qchist){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  for(auto hist : tthist_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : dbhist_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : gmhist_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vlhist_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : vllhist_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  for(auto hist : qchist_2D){if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

  if(doShapeSystematics){

    for(auto hist : dbhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vlhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vllhist_bUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_bDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}


    for(auto hist : dbhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vlhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vllhist_bUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_bDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  }


  outfile->cd();
  for(auto hist :  dthist) hist->Write();
  for(auto hist :  tthist) hist->Write();
  for(auto hist :  tthist_matched) hist->Write();
  for(auto hist :  tthist_unmatched) hist->Write();
  for(auto hist :  dbhist) hist->Write();
  for(auto hist :  gmhist) hist->Write();
  for(auto hist :  qchist) hist->Write();
  for(auto hist :  vlhist) hist->Write();
  for(auto hist :  vllhist) hist->Write();

  if(doShapeSystematics){
    for(auto hist :  dbhist_bUp) hist->Write();
    for(auto hist :  dbhist_bDw) hist->Write();
    for(auto hist :  dbhist_metJetUp) hist->Write();
    for(auto hist :  dbhist_metJetDw) hist->Write();
    for(auto hist :  dbhist_metResUp) hist->Write();
    for(auto hist :  dbhist_metResDw) hist->Write();
    for(auto hist :  dbhist_metUncUp) hist->Write();
    for(auto hist :  dbhist_metUncDw) hist->Write();

    for(auto hist :  gmhist_bUp) hist->Write();
    for(auto hist :  gmhist_bDw) hist->Write();
    for(auto hist :  gmhist_metJetUp) hist->Write();
    for(auto hist :  gmhist_metJetDw) hist->Write();
    for(auto hist :  gmhist_metResUp) hist->Write();
    for(auto hist :  gmhist_metResDw) hist->Write();
    for(auto hist :  gmhist_metUncUp) hist->Write();
    for(auto hist :  gmhist_metUncDw) hist->Write();

    for(auto hist :  vlhist_bUp) hist->Write();
    for(auto hist :  vlhist_bDw) hist->Write();
    for(auto hist :  vlhist_metJetUp) hist->Write();
    for(auto hist :  vlhist_metJetDw) hist->Write();
    for(auto hist :  vlhist_metResUp) hist->Write();
    for(auto hist :  vlhist_metResDw) hist->Write();
    for(auto hist :  vlhist_metUncUp) hist->Write();
    for(auto hist :  vlhist_metUncDw) hist->Write();

    for(auto hist :  vllhist_bUp) hist->Write();
    for(auto hist :  vllhist_bDw) hist->Write();
    for(auto hist :  vllhist_metJetUp) hist->Write();
    for(auto hist :  vllhist_metJetDw) hist->Write();
    for(auto hist :  vllhist_metResUp) hist->Write();
    for(auto hist :  vllhist_metResDw) hist->Write();
    for(auto hist :  vllhist_metUncUp) hist->Write();
    for(auto hist :  vllhist_metUncDw) hist->Write();
  }

  for(auto hist_2D : vlhist_2D) hist_2D->Write();
  for(auto hist_2D : vllhist_2D) hist_2D->Write();
  for(auto hist_2D : tthist_2D) hist_2D->Write();
  for(auto hist_2D : dbhist_2D) hist_2D->Write();
  for(auto hist_2D : gmhist_2D) hist_2D->Write();
  for(auto hist_2D : qchist_2D) hist_2D->Write();
  for(auto hist_2D : dthist_2D) hist_2D->Write();

  for(auto hist_2D : vlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : vllhist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : tthist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : dbhist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : gmhist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : qchist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : dthist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}

  if(doShapeSystematics){
    for(auto hist :  dbhist_bUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_bDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist :  gmhist_bUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_bDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist :  vlhist_bUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_bDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist :  vllhist_bUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_bDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
  }

  dtfile->Close();
  vlfile->Close();
  gmfile->Close();
  vllfile->Close();
  ttfile->Close();
  if(ttfile_alt) ttfile_alt->Close();
  dbfile->Close();
  qcfile->Close();

  dthist.clear();
  tthist.clear();
  tthist_alt.clear();
  qchist.clear();
  dbhist.clear();
  dbhist_bUp.clear();
  dbhist_bDw.clear();
  dbhist_metJetUp.clear();
  dbhist_metJetDw.clear();
  dbhist_metResUp.clear();
  dbhist_metResDw.clear();
  dbhist_metUncUp.clear();
  dbhist_metUncDw.clear();

  gmhist.clear();
  gmhist_bUp.clear();
  gmhist_bDw.clear();
  gmhist_metJetUp.clear();
  gmhist_metJetDw.clear();
  gmhist_metResUp.clear();
  gmhist_metResDw.clear();
  gmhist_metUncUp.clear();
  gmhist_metUncDw.clear();

  vlhist.clear();
  vlhist_bUp.clear();
  vlhist_bDw.clear();
  vlhist_metJetUp.clear();
  vlhist_metJetDw.clear();
  vlhist_metResUp.clear();
  vlhist_metResDw.clear();
  vlhist_metUncUp.clear();
  vlhist_metUncDw.clear();
  vllhist.clear();
  vllhist_bUp.clear();
  vllhist_bDw.clear();
  vllhist_metJetUp.clear();
  vllhist_metJetDw.clear();
  vllhist_metResUp.clear();
  vllhist_metResDw.clear();
  vllhist_metUncUp.clear();
  vllhist_metUncDw.clear();

  tthist_matched.clear();
  tthist_matched_alt.clear();
  tthist_unmatched.clear();
  tthist_unmatched_alt.clear();

  dthist_2D.clear();
  tthist_2D.clear();
  tthist_alt_2D.clear();
  qchist_2D.clear();
  dbhist_2D.clear();
  dbhist_bUp_2D.clear();
  dbhist_bDw_2D.clear();
  dbhist_metJetUp_2D.clear();
  dbhist_metJetDw_2D.clear();
  dbhist_metResUp_2D.clear();
  dbhist_metResDw_2D.clear();
  dbhist_metUncUp_2D.clear();
  dbhist_metUncDw_2D.clear();

  gmhist_2D.clear();
  gmhist_bUp_2D.clear();
  gmhist_bDw_2D.clear();
  gmhist_metJetUp_2D.clear();
  gmhist_metJetDw_2D.clear();
  gmhist_metResUp_2D.clear();
  gmhist_metResDw_2D.clear();
  gmhist_metUncUp_2D.clear();
  gmhist_metUncDw_2D.clear();

  vlhist_2D.clear();
  vlhist_bUp_2D.clear();
  vlhist_bDw_2D.clear();
  vlhist_metJetUp_2D.clear();
  vlhist_metJetDw_2D.clear();
  vlhist_metResUp_2D.clear();
  vlhist_metResDw_2D.clear();
  vlhist_metUncUp_2D.clear();
  vlhist_metUncDw_2D.clear();
  vllhist_2D.clear();
  vllhist_bUp_2D.clear();
  vllhist_bDw_2D.clear();
  vllhist_metJetUp_2D.clear();
  vllhist_metJetDw_2D.clear();
  vllhist_metResUp_2D.clear();
  vllhist_metResDw_2D.clear();
  vllhist_metUncUp_2D.clear();
  vllhist_metUncDw_2D.clear();
  
  cout << "Templates for the top control region computed ..." << endl;
}

