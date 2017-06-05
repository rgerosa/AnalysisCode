#include "makehist.h"
#include "TChain.h"

using namespace std;

static bool splitWPerFlavor = false;

void sigdatamchist(TFile* outfile,
                   const Category & category,
                   vector<string> observables,
		   vector<string> observables_2D,
                   const double & lumi,
		   const SamplesNLO & nloSamples,
		   const bool & doShapeSystematics  = false,
		   const bool & doAlternativeTop    = false,
                   const bool & blind               = false,
		   const bool & isHInv              = false,
		   const bool & applyPFWeight       = false,
		   const bool & useTheoriestKFactors = false) {


  // Files for Znunu,Wlnu, Zll, top, qcd , diboson, signal, data                                                                                                            
  TChain* zntree = new TChain("tree/tree");
  TChain* wltree = new TChain("tree/tree");
  TChain* zltree = new TChain("tree/tree");
  TChain* tttree = new TChain("tree/tree");
  TChain* tttree_alt = NULL;
  if(doAlternativeTop)
    tttree_alt = new TChain("tree/tree");
  TChain* ditree = new TChain("tree/tree");
  TChain* gmtree = new TChain("tree/tree");
  TChain* qcdtree = new TChain("tree/tree");
  TChain* ewkwtree = new TChain("tree/tree");
  TChain* ewkztree = new TChain("tree/tree");
  TChain* vgamtree = new TChain("tree/tree");
  TChain* zgtree = new TChain("tree/tree");
  TChain* wgtree = new TChain("tree/tree");
  TChain* dttree = new TChain("tree/tree");

  zntree->Add((baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/*root").c_str());
  wltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/*root").c_str());
  zltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/sigfilter/*root").c_str());
  tttree->Add((baseInputTreePath+"/TopAlternative/sigfilter/*root").c_str());
  if(doAlternativeTop and tttree_alt != NULL){
    tttree_alt->Add((baseInputTreePath+"/Top/sigfilter/*root").c_str());
  }
  qcdtree->Add((baseInputTreePath+"/QCD/sigfilter/*root").c_str());
  ditree->Add((baseInputTreePath+"/DiBoson/sigfilter/*root").c_str());
  gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/sigfilter/*root").c_str());
  ewkztree->Add((baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/*root").c_str());
  ewkwtree->Add((baseInputTreePath+"/WJetsEWK/sigfilter/*root").c_str());
  zgtree->Add((baseInputTreePath+"/ZnunuGJets/sigfilter/*root").c_str());
  wgtree->Add((baseInputTreePath+"/WGJets/sigfilter/*root").c_str());
  dttree->Add((baseInputTreePath+"/MET_reMiniAOD/sigfilter/*root").c_str());
  
  // make met histograms                                                                                                                                                        
  vector<TH1*> znhist;
  vector<TH1*> wlhist;

  vector<TH1*> wlhist_mu;
  vector<TH1*> wlhist_el;
  vector<TH1*> wlhist_ta;

  vector<TH1*> zlhist;
  vector<TH1*> zlhist_metJetUp;
  vector<TH1*> zlhist_metJetDw;
  vector<TH1*> zlhist_metResUp;
  vector<TH1*> zlhist_metResDw;
  vector<TH1*> zlhist_metUncUp;
  vector<TH1*> zlhist_metUncDw;

  vector<TH1*> tthist;
  vector<TH1*> tthist_metJetUp;
  vector<TH1*> tthist_metJetDw;
  vector<TH1*> tthist_metResUp;
  vector<TH1*> tthist_metResDw;
  vector<TH1*> tthist_metUncUp;
  vector<TH1*> tthist_metUncDw;

  vector<TH1*> tthist_alt;
  vector<TH1*> tthist_alt_metJetUp;
  vector<TH1*> tthist_alt_metJetDw;
  vector<TH1*> tthist_alt_metResUp;
  vector<TH1*> tthist_alt_metResDw;
  vector<TH1*> tthist_alt_metUncUp;
  vector<TH1*> tthist_alt_metUncDw;

  vector<TH1*> dihist;
  vector<TH1*> dihist_metJetUp;
  vector<TH1*> dihist_metJetDw;
  vector<TH1*> dihist_metResUp;
  vector<TH1*> dihist_metResDw;
  vector<TH1*> dihist_metUncUp;
  vector<TH1*> dihist_metUncDw;

  vector<TH1*> gmhist;
  vector<TH1*> gmhist_metJetUp;
  vector<TH1*> gmhist_metJetDw;
  vector<TH1*> gmhist_metResUp;
  vector<TH1*> gmhist_metResDw;
  vector<TH1*> gmhist_metUncUp;
  vector<TH1*> gmhist_metUncDw;

  vector<TH1*> qcdhist;
  vector<TH1*> ewkwhist;
  vector<TH1*> ewkwlhist_mu;
  vector<TH1*> ewkwlhist_el;
  vector<TH1*> ewkwlhist_ta;
  vector<TH1*> ewkzhist;
  vector<TH1*> vghist;

  vector<TH1*> dthist;

  vector<double> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;

    TH1F* znhist_temp = new TH1F(("zinvhist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* wlhist_temp = new TH1F(("wjethist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* zlhist_temp = new TH1F(("zjethist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* tthist_temp = new TH1F(("tbkghist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dihist_temp = new TH1F(("dbkghist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* gmhist_temp = new TH1F(("gbkghist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* qcdhist_temp = new TH1F(("qbkghist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkwhist_temp = new TH1F(("ewkbkgwhist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkzhist_temp = new TH1F(("ewkbkgzhist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vghist_temp = new TH1F(("vgbkghist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dthist_temp = new TH1F(("datahist_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

    znhist.push_back(dynamic_cast<TH1*>(znhist_temp));
    wlhist.push_back(dynamic_cast<TH1*>(wlhist_temp));
    if(splitWPerFlavor){
      TH1F* wlhist_mu_temp = new TH1F(("wjethist_mu_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* wlhist_el_temp = new TH1F(("wjethist_el_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* wlhist_ta_temp = new TH1F(("wjethist_ta_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      wlhist_mu.push_back(dynamic_cast<TH1*>(wlhist_mu_temp));
      wlhist_el.push_back(dynamic_cast<TH1*>(wlhist_el_temp));
      wlhist_ta.push_back(dynamic_cast<TH1*>(wlhist_ta_temp));

      TH1F* ewkwlhist_mu_temp = new TH1F(("ewkwjethist_mu_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* ewkwlhist_el_temp = new TH1F(("ewkwjethist_el_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* ewkwlhist_ta_temp = new TH1F(("ewkwjethist_ta_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      ewkwlhist_mu.push_back(dynamic_cast<TH1*>(ewkwlhist_mu_temp));
      ewkwlhist_el.push_back(dynamic_cast<TH1*>(ewkwlhist_el_temp));
      ewkwlhist_ta.push_back(dynamic_cast<TH1*>(ewkwlhist_ta_temp));
    }
    zlhist.push_back(dynamic_cast<TH1*>(zlhist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
    dihist.push_back(dynamic_cast<TH1*>(dihist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    ewkwhist.push_back(dynamic_cast<TH1*>(ewkwhist_temp));
    ewkzhist.push_back(dynamic_cast<TH1*>(ewkzhist_temp));
    vghist.push_back(dynamic_cast<TH1*>(vghist_temp));
    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));

    if(doAlternativeTop){
      TH1F* tthist_alt_temp = new TH1F(("tbkghist_alt_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      tthist_alt.push_back(dynamic_cast<TH1*>(tthist_alt_temp));
    }

    if(doShapeSystematics){


      TH1F* tthist_metJetUp_temp = new TH1F(("tbkghist_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metJetDw_temp = new TH1F(("tbkghist_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metResUp_temp = new TH1F(("tbkghist_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metResDw_temp = new TH1F(("tbkghist_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metUncUp_temp = new TH1F(("tbkghist_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metUncDw_temp = new TH1F(("tbkghist_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      tthist_metJetUp.push_back(dynamic_cast<TH1*>(tthist_metJetUp_temp));
      tthist_metJetDw.push_back(dynamic_cast<TH1*>(tthist_metJetDw_temp));
      tthist_metResUp.push_back(dynamic_cast<TH1*>(tthist_metResUp_temp));
      tthist_metResDw.push_back(dynamic_cast<TH1*>(tthist_metResDw_temp));
      tthist_metUncUp.push_back(dynamic_cast<TH1*>(tthist_metUncUp_temp));
      tthist_metUncDw.push_back(dynamic_cast<TH1*>(tthist_metUncDw_temp));

      TH1F* dihist_metJetUp_temp = new TH1F(("dbkghist_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dihist_metJetDw_temp = new TH1F(("dbkghist_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dihist_metResUp_temp = new TH1F(("dbkghist_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dihist_metResDw_temp = new TH1F(("dbkghist_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dihist_metUncUp_temp = new TH1F(("dbkghist_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dihist_metUncDw_temp = new TH1F(("dbkghist_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      dihist_metJetUp.push_back(dynamic_cast<TH1*>(dihist_metJetUp_temp));
      dihist_metJetDw.push_back(dynamic_cast<TH1*>(dihist_metJetDw_temp));
      dihist_metResUp.push_back(dynamic_cast<TH1*>(dihist_metResUp_temp));
      dihist_metResDw.push_back(dynamic_cast<TH1*>(dihist_metResDw_temp));
      dihist_metUncUp.push_back(dynamic_cast<TH1*>(dihist_metUncUp_temp));
      dihist_metUncDw.push_back(dynamic_cast<TH1*>(dihist_metUncDw_temp));

      TH1F* gmhist_metJetUp_temp = new TH1F(("gbkghist_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metJetDw_temp = new TH1F(("gbkghist_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metResUp_temp = new TH1F(("gbkghist_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metResDw_temp = new TH1F(("gbkghist_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metUncUp_temp = new TH1F(("gbkghist_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metUncDw_temp = new TH1F(("gbkghist_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      gmhist_metJetUp.push_back(dynamic_cast<TH1*>(gmhist_metJetUp_temp));
      gmhist_metJetDw.push_back(dynamic_cast<TH1*>(gmhist_metJetDw_temp));
      gmhist_metResUp.push_back(dynamic_cast<TH1*>(gmhist_metResUp_temp));
      gmhist_metResDw.push_back(dynamic_cast<TH1*>(gmhist_metResDw_temp));
      gmhist_metUncUp.push_back(dynamic_cast<TH1*>(gmhist_metUncUp_temp));
      gmhist_metUncDw.push_back(dynamic_cast<TH1*>(gmhist_metUncDw_temp));

      TH1F* zlhist_metJetUp_temp = new TH1F(("zjethist_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* zlhist_metJetDw_temp = new TH1F(("zjethist_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* zlhist_metResUp_temp = new TH1F(("zjethist_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* zlhist_metResDw_temp = new TH1F(("zjethist_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* zlhist_metUncUp_temp = new TH1F(("zjethist_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* zlhist_metUncDw_temp = new TH1F(("zjethist_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      zlhist_metJetUp.push_back(dynamic_cast<TH1*>(zlhist_metJetUp_temp));
      zlhist_metJetDw.push_back(dynamic_cast<TH1*>(zlhist_metJetDw_temp));
      zlhist_metResUp.push_back(dynamic_cast<TH1*>(zlhist_metResUp_temp));
      zlhist_metResDw.push_back(dynamic_cast<TH1*>(zlhist_metResDw_temp));
      zlhist_metUncUp.push_back(dynamic_cast<TH1*>(zlhist_metUncUp_temp));
      zlhist_metUncDw.push_back(dynamic_cast<TH1*>(zlhist_metUncDw_temp));
      
      if(doAlternativeTop){


	TH1F* tthist_alt_metJetUp_temp = new TH1F(("tbkghist_alt_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* tthist_alt_metJetDw_temp = new TH1F(("tbkghist_alt_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* tthist_alt_metResUp_temp = new TH1F(("tbkghist_alt_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* tthist_alt_metResDw_temp = new TH1F(("tbkghist_alt_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* tthist_alt_metUncUp_temp = new TH1F(("tbkghist_alt_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	TH1F* tthist_alt_metUncDw_temp = new TH1F(("tbkghist_alt_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
	
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
  vector<TH2*> wlhist_mu_2D;
  vector<TH2*> wlhist_el_2D;
  vector<TH2*> wlhist_ta_2D;
  vector<TH2*> zlhist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_alt_2D;
  vector<TH2*> dihist_2D;
  vector<TH2*> gmhist_2D;
  vector<TH2*> vghist_2D;
  vector<TH2*> qcdhist_2D;

  vector<TH2*> zlhist_metJetUp_2D;
  vector<TH2*> zlhist_metJetDw_2D;
  vector<TH2*> zlhist_metResUp_2D;
  vector<TH2*> zlhist_metResDw_2D;
  vector<TH2*> zlhist_metUncUp_2D;
  vector<TH2*> zlhist_metUncDw_2D;
  vector<TH2*> tthist_metJetUp_2D;
  vector<TH2*> tthist_metJetDw_2D;
  vector<TH2*> tthist_metResUp_2D;
  vector<TH2*> tthist_metResDw_2D;
  vector<TH2*> tthist_metUncUp_2D;
  vector<TH2*> tthist_metUncDw_2D;
  vector<TH2*> tthist_alt_metJetUp_2D;
  vector<TH2*> tthist_alt_metJetDw_2D;
  vector<TH2*> tthist_alt_metResUp_2D;
  vector<TH2*> tthist_alt_metResDw_2D;
  vector<TH2*> tthist_alt_metUncUp_2D;
  vector<TH2*> tthist_alt_metUncDw_2D;
  vector<TH2*> dihist_metJetUp_2D;
  vector<TH2*> dihist_metJetDw_2D;
  vector<TH2*> dihist_metResUp_2D;
  vector<TH2*> dihist_metResDw_2D;
  vector<TH2*> dihist_metUncUp_2D;
  vector<TH2*> dihist_metUncDw_2D;
  vector<TH2*> gmhist_metJetUp_2D;
  vector<TH2*> gmhist_metJetDw_2D;
  vector<TH2*> gmhist_metResUp_2D;
  vector<TH2*> gmhist_metResDw_2D;
  vector<TH2*> gmhist_metUncUp_2D;
  vector<TH2*> gmhist_metUncDw_2D;
  vector<TH2*> ewkwhist_2D;
  vector<TH2*> ewkwlhist_mu_2D;
  vector<TH2*> ewkwlhist_el_2D;
  vector<TH2*> ewkwlhist_ta_2D;
  vector<TH2*> ewkzhist_2D;

  vector<TH2*> dthist_2D;


  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable "+obs+"--> please define it"<<endl;
    
    TH2F* znhist_temp = new TH2F(("zinvhist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* wlhist_temp = new TH2F(("wjethist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* zlhist_temp = new TH2F(("zjethist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* tthist_temp = new TH2F(("tbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dihist_temp = new TH2F(("dbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* gmhist_temp = new TH2F(("gbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* qcdhist_temp = new TH2F(("qbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkwhist_temp = new TH2F(("ewkwbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkzhist_temp = new TH2F(("ewkzbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vghist_temp = new TH2F(("vgbkghist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dthist_temp = new TH2F(("datahist_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

    znhist_2D.push_back(dynamic_cast<TH2*>(znhist_temp));
    wlhist_2D.push_back(dynamic_cast<TH2*>(wlhist_temp));
    if(splitWPerFlavor){
      TH2F* wlhist_mu_temp = new TH2F(("wjethist_mu_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* wlhist_el_temp = new TH2F(("wjethist_el_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* wlhist_ta_temp = new TH2F(("wjethist_ta_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      wlhist_mu_2D.push_back(dynamic_cast<TH2*>(wlhist_mu_temp));
      wlhist_el_2D.push_back(dynamic_cast<TH2*>(wlhist_el_temp));
      wlhist_ta_2D.push_back(dynamic_cast<TH2*>(wlhist_ta_temp));

      TH2F* ewkwlhist_mu_temp = new TH2F(("ewkwjethist_mu_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* ewkwlhist_el_temp = new TH2F(("ewkwjethist_el_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* ewkwlhist_ta_temp = new TH2F(("ewkwjethist_ta_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      ewkwlhist_mu_2D.push_back(dynamic_cast<TH2*>(ewkwlhist_mu_temp));
      ewkwlhist_el_2D.push_back(dynamic_cast<TH2*>(ewkwlhist_el_temp));
      ewkwlhist_ta_2D.push_back(dynamic_cast<TH2*>(ewkwlhist_ta_temp));
    }
    zlhist_2D.push_back(dynamic_cast<TH2*>(zlhist_temp));
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    qcdhist_2D.push_back(dynamic_cast<TH2*>(qcdhist_temp));
    dihist_2D.push_back(dynamic_cast<TH2*>(dihist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    ewkwhist_2D.push_back(dynamic_cast<TH2*>(ewkwhist_temp));
    ewkzhist_2D.push_back(dynamic_cast<TH2*>(ewkzhist_temp));
    vghist_2D.push_back(dynamic_cast<TH2*>(vghist_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));

    if(doAlternativeTop){
      TH2F* tthist_alt_temp = new TH2F(("tbkghist_alt_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      tthist_alt_2D.push_back(dynamic_cast<TH2*>(tthist_alt_temp));
    }    


    if(doShapeSystematics){

      TH2F* zlhist_metJetUp_temp = new TH2F(("zjethist_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* zlhist_metJetDw_temp = new TH2F(("zjethist_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* zlhist_metResUp_temp = new TH2F(("zjethist_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* zlhist_metResDw_temp = new TH2F(("zjethist_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* zlhist_metUncUp_temp = new TH2F(("zjethist_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* zlhist_metUncDw_temp = new TH2F(("zjethist_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      zlhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(zlhist_metJetUp_temp));
      zlhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(zlhist_metJetDw_temp));
      zlhist_metResUp_2D.push_back(dynamic_cast<TH2*>(zlhist_metResUp_temp));
      zlhist_metResDw_2D.push_back(dynamic_cast<TH2*>(zlhist_metResDw_temp));
      zlhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(zlhist_metUncUp_temp));
      zlhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(zlhist_metUncDw_temp));

      TH2F* tthist_metJetUp_temp = new TH2F(("tbkghist_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metJetDw_temp = new TH2F(("tbkghist_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metResUp_temp = new TH2F(("tbkghist_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metResDw_temp = new TH2F(("tbkghist_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metUncUp_temp = new TH2F(("tbkghist_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metUncDw_temp = new TH2F(("tbkghist_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      tthist_metJetUp_2D.push_back(dynamic_cast<TH2*>(tthist_metJetUp_temp));
      tthist_metJetDw_2D.push_back(dynamic_cast<TH2*>(tthist_metJetDw_temp));
      tthist_metResUp_2D.push_back(dynamic_cast<TH2*>(tthist_metResUp_temp));
      tthist_metResDw_2D.push_back(dynamic_cast<TH2*>(tthist_metResDw_temp));
      tthist_metUncUp_2D.push_back(dynamic_cast<TH2*>(tthist_metUncUp_temp));
      tthist_metUncDw_2D.push_back(dynamic_cast<TH2*>(tthist_metUncDw_temp));

      if(doAlternativeTop){

	TH2F* tthist_alt_metJetUp_temp = new TH2F(("tbkghist_alt_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* tthist_alt_metJetDw_temp = new TH2F(("tbkghist_alt_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* tthist_alt_metResUp_temp = new TH2F(("tbkghist_alt_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* tthist_alt_metResDw_temp = new TH2F(("tbkghist_alt_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* tthist_alt_metUncUp_temp = new TH2F(("tbkghist_alt_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	TH2F* tthist_alt_metUncDw_temp = new TH2F(("tbkghist_alt_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
	tthist_alt_metJetUp_2D.push_back(dynamic_cast<TH2*>(tthist_alt_metJetUp_temp));
	tthist_alt_metJetDw_2D.push_back(dynamic_cast<TH2*>(tthist_alt_metJetDw_temp));
	tthist_alt_metResUp_2D.push_back(dynamic_cast<TH2*>(tthist_alt_metResUp_temp));
	tthist_alt_metResDw_2D.push_back(dynamic_cast<TH2*>(tthist_alt_metResDw_temp));
	tthist_alt_metUncUp_2D.push_back(dynamic_cast<TH2*>(tthist_alt_metUncUp_temp));
	tthist_alt_metUncDw_2D.push_back(dynamic_cast<TH2*>(tthist_alt_metUncDw_temp));

      }

      TH2F* dihist_metJetUp_temp = new TH2F(("dbkghist_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dihist_metJetDw_temp = new TH2F(("dbkghist_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dihist_metResUp_temp = new TH2F(("dbkghist_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dihist_metResDw_temp = new TH2F(("dbkghist_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dihist_metUncUp_temp = new TH2F(("dbkghist_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dihist_metUncDw_temp = new TH2F(("dbkghist_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      dihist_metJetUp_2D.push_back(dynamic_cast<TH2*>(dihist_metJetUp_temp));
      dihist_metJetDw_2D.push_back(dynamic_cast<TH2*>(dihist_metJetDw_temp));
      dihist_metResUp_2D.push_back(dynamic_cast<TH2*>(dihist_metResUp_temp));
      dihist_metResDw_2D.push_back(dynamic_cast<TH2*>(dihist_metResDw_temp));
      dihist_metUncUp_2D.push_back(dynamic_cast<TH2*>(dihist_metUncUp_temp));
      dihist_metUncDw_2D.push_back(dynamic_cast<TH2*>(dihist_metUncDw_temp));

      TH2F* gmhist_metJetUp_temp = new TH2F(("gbkghist_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metJetDw_temp = new TH2F(("gbkghist_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metResUp_temp = new TH2F(("gbkghist_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metResDw_temp = new TH2F(("gbkghist_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metUncUp_temp = new TH2F(("gbkghist_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metUncDw_temp = new TH2F(("gbkghist_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      gmhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetUp_temp));
      gmhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetDw_temp));
      gmhist_metResUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metResUp_temp));
      gmhist_metResDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metResDw_temp));
      gmhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncUp_temp));
      gmhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncDw_temp));

    }    
  }

  
  // get k-factors NLO                                                                                                                                                           
  TFile* kffile = NULL;
  TFile* kffile_alt = NULL;
  TFile* kffile_alt2 = NULL;
  TFile* kffile_unlops = NULL;

  TH1*   aunlopshist = NULL;
  TH1*  znlohist = NULL;
  TH1*  zlohist  = NULL;
  TH1* zewkhist  = NULL;
  TH1*  wnlohist = NULL;
  TH1*  wlohist  = NULL;
  TH1* wewkhist  = NULL;
  TH1* anlohist  = NULL;
  TH1*  alohist  = NULL;
  TH1* aewkhist  = NULL;
  TH1* reweight_zvv =  NULL;
  TH1* reweight_zll =  NULL;
  TH1* reweight_wln =  NULL;

  if(not useTheoriestKFactors){
    
    kffile = TFile::Open(kfactorFile.c_str());
    znlohist = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
    zlohist  = (TH1*) kffile->Get("ZJets_LO/inv_pt");
    zewkhist  = (TH1*) kffile->Get("EWKcorr/Z");
    
    if(zewkhist)
      zewkhist->Divide(znlohist);
    if(znlohist)
      znlohist->Divide(zlohist);

    wnlohist = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
    wlohist  = (TH1*) kffile->Get("WJets_LO/inv_pt");
    wewkhist  = (TH1*) kffile->Get("EWKcorr/W");
    
    if(wewkhist)
      wewkhist->Divide(wnlohist);
    if(wnlohist)
      wnlohist->Divide(wlohist);
    
    kffile_unlops = TFile::Open(kfactorFileUNLOPS.c_str());
    aunlopshist = (TH1*) kffile_unlops->Get("Func");
    anlohist  = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    alohist  = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    aewkhist  = (TH1*) kffile->Get("EWKcorr/photon");
    
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);
  }
  else{


    kffile = TFile::Open(kFactorTheoristFile_zvv.c_str());
    TH1* kfact_nloqcd_zvv = (TH1*) kffile->Get("vvj_pTV_K_NLO");
    TH1* kfact_nloewk_zvv = (TH1*) kffile->Get("vvj_pTV_kappa_NLO_EW");
    TH1* kfact_sudewk_zvv = (TH1*) kffile->Get("vvj_pTV_kappa_NNLO_Sud");
    reweight_zvv = (TH1*) kfact_nloqcd_zvv->Clone("reweight_zvv");
    reweight_zvv->Reset("ICES");

    for(int iBin = 1; iBin <= reweight_zvv->GetNbinsX(); iBin++)
      reweight_zvv->SetBinContent(iBin,kfact_nloqcd_zvv->GetBinContent(iBin)*(1+kfact_nloewk_zvv->GetBinContent(iBin)+kfact_sudewk_zvv->GetBinContent(iBin)));

    kffile_alt = TFile::Open(kFactorTheoristFile_zll.c_str());
    TH1* kfact_nloqcd_zll = (TH1*) kffile_alt->Get("eej_pTV_K_NLO");
    TH1* kfact_nloewk_zll = (TH1*) kffile_alt->Get("eej_pTV_kappa_NLO_EW");
    TH1* kfact_sudewk_zll = (TH1*) kffile_alt->Get("eej_pTV_kappa_NNLO_Sud");
    reweight_zll = (TH1*) kfact_nloqcd_zll->Clone("reweight_zll");
    reweight_zll->Reset("ICES");

    for(int iBin = 1; iBin < reweight_zll->GetNbinsX(); iBin++)
      reweight_zll->SetBinContent(iBin,kfact_nloqcd_zll->GetBinContent(iBin)*(1+kfact_nloewk_zll->GetBinContent(iBin)+kfact_sudewk_zll->GetBinContent(iBin)));

    kffile_alt2 = TFile::Open(kFactorTheoristFile_wln.c_str());
    TH1* kfact_nloqcd_wln = (TH1*) kffile_alt2->Get("evj_pTV_K_NLO");
    TH1* kfact_nloewk_wln = (TH1*) kffile_alt2->Get("evj_pTV_kappa_NLO_EW");
    TH1* kfact_sudewk_wln = (TH1*) kffile_alt2->Get("evj_pTV_kappa_NNLO_Sud");
    reweight_wln = (TH1*) kfact_nloqcd_wln->Clone("reweight_wln");
    reweight_wln->Reset("ICES");

    for(int iBin = 1; iBin < reweight_wln->GetNbinsX(); iBin++)
      reweight_wln->SetBinContent(iBin,kfact_nloqcd_wln->GetBinContent(iBin)*(1+kfact_nloewk_wln->GetBinContent(iBin)+kfact_sudewk_wln->GetBinContent(iBin)));
    
    // no gamma+jets k-factors yet
    anlohist  = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    alohist  = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    aewkhist  = (TH1*) kffile->Get("EWKcorr/photon");
    
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);
    
  }
  
  // apply EWK and QCD corrections                                                                                                                                              
  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> dyhists;
  vector<TH1*> whists;
  vector<TH1*> ahists;

  if(not useTheoriestKFactors){
    zhists.push_back(znlohist); 
    zhists.push_back(zewkhist);
    dyhists.push_back(znlohist); 
    dyhists.push_back(zewkhist);
    whists.push_back(wnlohist); 
    whists.push_back(wewkhist);    
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);
    ahists.push_back(aunlopshist);
  }
  else{
    zhists.push_back(reweight_zvv);
    dyhists.push_back(reweight_zll);
    whists.push_back(reweight_wln);
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);
  }

  float scale = 1;
  if(nloSamples.useWJetsNLO){
    // temp fix since we always use the W+jets NLO
    whists.clear();
    whists.push_back(wewkhist);
  }
  if(nloSamples.useZJetsNLO){
    zhists.clear();
    zhists.push_back(zewkhist);
    scale = 3.;
  }
  if(nloSamples.useDYJetsNLO){
    dyhists.clear();
    dyhists.push_back(zewkhist);
  }
  if(nloSamples.usePhotonJetsNLO){
    ahists.clear();
    ahists.push_back(aewkhist);
    ahists.push_back(aunlopshist);
  }

  // special corrections for VBF topology
  TFile* kfactzjet_vbf = NULL;
  TFile* kfactwjet_vbf = NULL;
  TFile* kfactgjet_vbf = NULL;

  if((category == Category::VBF or category == Category::VBFrelaxed) and not useTheoriestKFactors){ // apply further k-factors going to the VBF selections

    /////////////
    kfactzjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_zjets_v2.root");
    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets_v2.root");
    kfactgjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_gjets_v2.root");
    
    ////////////
    TH1* zjet_nlo_vbf = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      zjet_nlo_vbf = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_vbf_relaxed");
    TH1* zjet_nlo_mj  = (TH1*) kfactzjet_vbf->Get("bosonPt_NLO_monojet");

    if(category == Category::VBF)
      zjet_nlo_vbf->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      zjet_nlo_vbf->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
    zjet_nlo_mj->Divide((TH1*) kfactzjet_vbf->Get("bosonPt_LO_monojet"));
    zjet_nlo_vbf->Divide(zjet_nlo_mj);

    if(not nloSamples.useZJetsNLO)
      zhists.push_back(zjet_nlo_vbf);
    if(not nloSamples.useDYJetsNLO)
      dyhists.push_back(zjet_nlo_vbf);
    
    ///////////////////
    TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf_relaxed");
    TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");

    if(category == Category::VBF)
      wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
    
    wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);
    if(not nloSamples.useWJetsNLO)
      whists.push_back(wjet_nlo_vbf);

    TH1* gjet_nlo_vbf = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      gjet_nlo_vbf = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_vbf_relaxed");

    TH1* gjet_nlo_mj  = (TH1*) kfactgjet_vbf->Get("bosonPt_NLO_monojet");
    if(category == Category::VBF)
      gjet_nlo_vbf->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      gjet_nlo_vbf->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
      
    gjet_nlo_mj->Divide((TH1*) kfactgjet_vbf->Get("bosonPt_LO_monojet"));
    gjet_nlo_vbf->Divide(gjet_nlo_mj);

    if(not nloSamples.usePhotonJetsNLO)
      ahists.push_back(gjet_nlo_vbf);
  }

  // for electroweak backgrounds
  TFile* kffile_zewk = TFile::Open(kFactorFile_zjetewk.c_str(),"READ");
  TFile* kffile_wewk = TFile::Open(kFactorFile_wjetewk.c_str(),"READ");
  vector<TH2*> zewkhists;
  vector<TH2*> wewkhists;
  vector<TH2*> ehists2D;
  if(applyEWKVKfactor){
    zewkhists.push_back((TH2*) kffile_zewk->Get("TH2F_kFactor"));
    wewkhists.push_back((TH2*) kffile_wewk->Get("TH2F_kFactor"));
  }

  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;

  cout<<"signal region: Z->nunu sample "<<endl;

  if(not nloSamples.useZJetsNLO)
    makehist4(zntree,znhist,znhist_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  else
    makehist4(zntree,znhist,znhist_2D,true,Sample::sig,category,false,scale,lumi,zhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  cout<<"signal region: W+jets sample "<<endl;
  makehist4(wltree,wlhist,wlhist_2D,true,Sample::sig,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  if(splitWPerFlavor){
    cout<<"signal region: W+jets sample: wmn-gen "<<endl;
    makehist4(wltree,wlhist_mu,wlhist_mu_2D,true,Sample::sig,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight,ehists2D,-1,NULL,NULL,NULL,"muon");
    cout<<"signal region: W+jets sample  wen-gen "<<endl;
    makehist4(wltree,wlhist_el,wlhist_el_2D,true,Sample::sig,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight,ehists2D,-1,NULL,NULL,NULL,"electron");
    cout<<"signal region: W+jets sample  wta-gen "<<endl;
    makehist4(wltree,wlhist_ta,wlhist_ta_2D,true,Sample::sig,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight,ehists2D,-1,NULL,NULL,NULL,"tau");
  }
  cout<<"signal region: Z+jets sample "<<endl;
  makehist4(zltree,zlhist,zlhist_2D,true,Sample::sig,category,false,1.00,lumi,dyhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"signal region: gamma+jets sample "<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"signal region: ewkw+jets sample "<<endl;
  makehist4(ewkwtree,ewkwhist,ewkwhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight,wewkhists);
  if(splitWPerFlavor){
    cout<<"signal region: ewkw+jets sample: wmn-gen "<<endl;
    makehist4(ewkwtree,ewkwlhist_mu,ewkwlhist_mu_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight,wewkhists,-1,NULL,NULL,NULL,"muon");
    cout<<"signal region: ewkw+jets sample  wen-gen "<<endl;
    makehist4(ewkwtree,ewkwlhist_el,ewkwlhist_el_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight,wewkhists,-1,NULL,NULL,NULL,"electron");
    cout<<"signal region: ewkw+jets sample  wta-gen "<<endl;
    makehist4(ewkwtree,ewkwlhist_ta,ewkwlhist_ta_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight,wewkhists,-1,NULL,NULL,NULL,"tau");
  }
  cout<<"signal region: ewkz+jets sample "<<endl;
  makehist4(ewkztree,ewkzhist,ewkzhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight,zewkhists); // temp fix for a wrong xsec
  cout<<"signal region: TTbar sample "<<endl;
  makehist4(tttree,tthist,tthist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    //alternative ttbar             
  if(doAlternativeTop){
    cout<<"signal region: TTbar alternative sample "<<endl;
    makehist4(tttree_alt,tthist_alt,tthist_alt_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }

  cout<<"signal region: Diboson sample "<<endl;  
  makehist4(ditree,dihist,dihist_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"signal region: Wgamma+jets"<<endl;
  makehist4(wgtree,vghist,vghist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"signal region: Zgamma+jets"<<endl;
  makehist4(zgtree,vghist,vghist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  cout<<"signal region: QCD sample "<<endl;
  makehist4(qcdtree,qcdhist,qcdhist_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  if(doShapeSystematics){
    cout<<"signal region analysis --> do top shape sys "<<endl;
    makehist4(tttree,tthist_metJetUp,tthist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metResUp,tthist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metResDw,tthist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metUncUp,tthist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metUncDw,tthist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    
    if(doAlternativeTop){
      cout<<"signal region analysis --> do top alternative shape sys "<< endl;
      makehist4(tttree_alt,tthist_alt_metJetUp,tthist_alt_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
      makehist4(tttree_alt,tthist_alt_metJetDw,tthist_alt_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
      makehist4(tttree_alt,tthist_alt_metResUp,tthist_alt_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
      makehist4(tttree_alt,tthist_alt_metResDw,tthist_alt_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
      makehist4(tttree_alt,tthist_alt_metUncUp,tthist_alt_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
      makehist4(tttree_alt,tthist_alt_metUncDw,tthist_alt_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    }


    cout<<"signal region analysis --> do diboson shape sys "<<endl;
    makehist4(ditree,dihist_metJetUp,dihist_metJetUp_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(ditree,dihist_metJetDw,dihist_metJetDw_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(ditree,dihist_metResUp,dihist_metResUp_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(ditree,dihist_metResDw,dihist_metResDw_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(ditree,dihist_metUncUp,dihist_metUncUp_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(ditree,dihist_metUncDw,dihist_metUncDw_2D,true,Sample::sig,category,isWJet,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"signal region analysis --> do DYJets shape sys "<<endl;
    makehist4(zltree,zlhist_metJetUp,zlhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(zltree,zlhist_metJetDw,zlhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(zltree,zlhist_metResUp,zlhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(zltree,zlhist_metResDw,zlhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(zltree,zlhist_metUncUp,zlhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(zltree,zlhist_metUncDw,zlhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,zhists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"signal region analysis --> do gamma+jets shape sys "<<endl;
    makehist4(gmtree,gmhist_metJetUp,gmhist_metJetUp_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metJetDw,gmhist_metJetDw_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metResUp,gmhist_metResUp_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metResDw,gmhist_metResDw_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metUncUp,gmhist_metUncUp_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metUncDw,gmhist_metUncDw_2D,true,Sample::sig,category,false,1.00,lumi,ahists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    
  }

  // take average of ttbar                                                                                                                                                     
  if(doAlternativeTop){

    for(size_t iHisto = 0; iHisto < tthist.size(); iHisto++)
      makeAverage(tthist.at(iHisto),tthist_alt.at(iHisto));      

    for(size_t iHisto = 0; iHisto < tthist_2D.size(); iHisto++)
      makeAverage(tthist_2D.at(iHisto),tthist_alt_2D.at(iHisto));

    if(doShapeSystematics){

      for(size_t iHisto = 0; iHisto < tthist_metJetUp.size(); iHisto++)
	makeAverage(tthist_metJetUp.at(iHisto),tthist_alt_metJetUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metJetUp_2D.size(); iHisto++)
	makeAverage(tthist_metJetUp_2D.at(iHisto),tthist_alt_metJetUp_2D.at(iHisto));

      for(size_t iHisto = 0; iHisto < tthist_metJetDw.size(); iHisto++)
	makeAverage(tthist_metJetDw.at(iHisto),tthist_alt_metJetDw.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metJetDw_2D.size(); iHisto++)
	makeAverage(tthist_metJetDw_2D.at(iHisto),tthist_alt_metJetDw_2D.at(iHisto));

      for(size_t iHisto = 0; iHisto < tthist_metResUp.size(); iHisto++)
	makeAverage(tthist_metResUp.at(iHisto),tthist_alt_metResUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metResUp_2D.size(); iHisto++)
	makeAverage(tthist_metResUp_2D.at(iHisto),tthist_alt_metResUp_2D.at(iHisto));

      for(size_t iHisto = 0; iHisto < tthist_metResDw.size(); iHisto++)
	makeAverage(tthist_metResDw.at(iHisto),tthist_alt_metResDw.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metResDw_2D.size(); iHisto++)
	makeAverage(tthist_metResDw_2D.at(iHisto),tthist_alt_metResDw_2D.at(iHisto));

      for(size_t iHisto = 0; iHisto < tthist_metUncUp.size(); iHisto++)
	makeAverage(tthist_metUncUp.at(iHisto),tthist_alt_metUncUp.at(iHisto));
      for(size_t iHisto = 0; iHisto < tthist_metUncUp_2D.size(); iHisto++)
	makeAverage(tthist_metUncUp_2D.at(iHisto),tthist_alt_metUncUp_2D.at(iHisto));

      for(size_t iHisto = 0; iHisto < tthist_metUncDw.size(); iHisto++)
	makeAverage(tthist_metUncDw.at(iHisto),tthist_alt_metUncDw.at(iHisto));	
      for(size_t iHisto = 0; iHisto < tthist_metUncDw_2D.size(); iHisto++)
	makeAverage(tthist_metUncDw_2D.at(iHisto),tthist_alt_metUncDw_2D.at(iHisto));	
    }
  }

  //smooth
  if(doSmoothing){
    for(auto hist : tthist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qcdhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    
    for(auto hist : tthist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dihist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : zlhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    //  for(auto hist : qcdhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    
    if(doShapeSystematics){
      
      for(auto hist : tthist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : dihist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : zlhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : gmhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : tthist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : dihist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dihist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : zlhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : zlhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : gmhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
    }
  }

  // data                                                
  cout<<"signal region analysis --> loop on data "<<endl;
  makehist4(dttree,dthist,dthist_2D,false,Sample::sig,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  
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
        dthist.at(ihist)->SetBinContent(i,int(binval));
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
  if(not outfile->GetDirectory("SR"))
    outfile->mkdir("SR");
  outfile->cd("SR");
  // store histograms                                                                                                                                                        
  for(auto hist : znhist) hist->Write();
  for(auto hist : wlhist) hist->Write();
  if(splitWPerFlavor){
    for(auto hist : wlhist_mu) hist->Write();
    for(auto hist : wlhist_el) hist->Write();
    for(auto hist : wlhist_ta) hist->Write();
  }
  for(auto hist : zlhist) hist->Write();
  for(auto hist : tthist) hist->Write();
  for(auto hist : dihist) hist->Write();
  for(auto hist : gmhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();
  for(auto hist : ewkwhist) hist->Write();
  if(splitWPerFlavor){
    for(auto hist : ewkwlhist_mu) hist->Write();
    for(auto hist : ewkwlhist_el) hist->Write();
    for(auto hist : ewkwlhist_ta) hist->Write();
  }
  for(auto hist : ewkzhist) hist->Write();
  for(auto hist : vghist) hist->Write();
  for(auto hist : dthist) hist->Write();

  //
  if(doShapeSystematics){
    
    outfile->cd();
    if(not outfile->GetDirectory("SR/sysShape"))
      outfile->mkdir("SR/sysShape");
    outfile->cd("SR/sysShape");
    for(auto hist : zlhist_metJetUp) hist->Write();
    for(auto hist : zlhist_metJetDw) hist->Write();
    for(auto hist : zlhist_metResUp) hist->Write();
    for(auto hist : zlhist_metResDw) hist->Write();
    for(auto hist : zlhist_metUncUp) hist->Write();
    for(auto hist : zlhist_metUncDw) hist->Write();

    for(auto hist : gmhist_metJetUp) hist->Write();
    for(auto hist : gmhist_metJetDw) hist->Write();
    for(auto hist : gmhist_metResUp) hist->Write();
    for(auto hist : gmhist_metResDw) hist->Write();
    for(auto hist : gmhist_metUncUp) hist->Write();
    for(auto hist : gmhist_metUncDw) hist->Write();

    for(auto hist : tthist_metJetUp) hist->Write();
    for(auto hist : tthist_metJetDw) hist->Write();
    for(auto hist : tthist_metResUp) hist->Write();
    for(auto hist : tthist_metResDw) hist->Write();
    for(auto hist : tthist_metUncUp) hist->Write();
    for(auto hist : tthist_metUncDw) hist->Write();

    for(auto hist : dihist_metJetUp) hist->Write();
    for(auto hist : dihist_metJetDw) hist->Write();
    for(auto hist : dihist_metResUp) hist->Write();
    for(auto hist : dihist_metResDw) hist->Write();
    for(auto hist : dihist_metUncUp) hist->Write();
    for(auto hist : dihist_metUncDw) hist->Write();
    
  }

  // store hist_2Dograms                                                                                                                                                     
  outfile->cd();
  outfile->cd("SR");
  for(auto hist_2D : znhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : wlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  if(splitWPerFlavor){
    for(auto hist_2D : wlhist_mu_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : wlhist_el_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : wlhist_ta_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  }
  for(auto hist_2D : zlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : gmhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : tthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dihist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : qcdhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkwhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  if(splitWPerFlavor){
    for(auto hist_2D : ewkwlhist_mu_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : ewkwlhist_el_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : ewkwlhist_ta_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  }
  for(auto hist_2D : ewkzhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vghist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }


  if(doShapeSystematics){
    
    outfile->cd();
    outfile->cd("SR/sysShape");
    for(auto hist_2D : zlhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : zlhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    for(auto hist_2D : gmhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : gmhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    for(auto hist_2D : tthist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : tthist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    for(auto hist_2D : dihist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
    for(auto hist_2D : dihist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

    
  }

  outfile->cd();

  zhists.clear();
  whists.clear();
  dyhists.clear();
  ahists.clear();
  zewkhists.clear();
  wewkhists.clear();

  for(auto hist: znhist) delete hist;
  for(auto hist: wlhist) delete hist;
  for(auto hist: wlhist_mu) delete hist;
  for(auto hist: wlhist_el) delete hist;
  for(auto hist: wlhist_ta) delete hist;
  for(auto hist: zlhist) delete hist;
  for(auto hist: zlhist_metJetUp) delete hist;
  for(auto hist: zlhist_metJetDw) delete hist;
  for(auto hist: zlhist_metResUp) delete hist;
  for(auto hist: zlhist_metResDw) delete hist;
  for(auto hist: zlhist_metUncUp) delete hist;
  for(auto hist: zlhist_metUncDw) delete hist;
  for(auto hist: tthist) delete hist;
  for(auto hist: tthist_metJetUp) delete hist;
  for(auto hist: tthist_metJetDw) delete hist;
  for(auto hist: tthist_metResUp) delete hist;
  for(auto hist: tthist_metResDw) delete hist;
  for(auto hist: tthist_metUncUp) delete hist;
  for(auto hist: tthist_metUncDw) delete hist;
  for(auto hist: tthist_alt) delete hist;
  for(auto hist: tthist_alt_metJetUp) delete hist;
  for(auto hist: tthist_alt_metJetDw) delete hist;
  for(auto hist: tthist_alt_metResUp) delete hist;
  for(auto hist: tthist_alt_metResDw) delete hist;
  for(auto hist: tthist_alt_metUncUp) delete hist;
  for(auto hist: tthist_alt_metUncDw) delete hist;
  for(auto hist: dihist) delete hist;
  for(auto hist: dihist_metJetUp) delete hist;
  for(auto hist: dihist_metJetDw) delete hist;
  for(auto hist: dihist_metResUp) delete hist;
  for(auto hist: dihist_metResDw) delete hist;
  for(auto hist: dihist_metUncUp) delete hist;
  for(auto hist: dihist_metUncDw) delete hist;
  for(auto hist: gmhist) delete hist;
  for(auto hist: gmhist_metJetUp) delete hist;
  for(auto hist: gmhist_metJetDw) delete hist;
  for(auto hist: gmhist_metResUp) delete hist;
  for(auto hist: gmhist_metResDw) delete hist;
  for(auto hist: gmhist_metUncUp) delete hist;
  for(auto hist: gmhist_metUncDw) delete hist;
  for(auto hist: qcdhist) delete hist;
  for(auto hist: ewkwhist) delete hist;
  for(auto hist: ewkwlhist_mu) delete hist;
  for(auto hist: ewkwlhist_el) delete hist;
  for(auto hist: ewkwlhist_ta) delete hist;
  for(auto hist: ewkzhist) delete hist;
  for(auto hist: vghist) delete hist;
  for(auto hist: dthist) delete hist;

  znhist.clear();
  wlhist.clear();
  wlhist_mu.clear();
  wlhist_el.clear();
  wlhist_ta.clear();
  zlhist.clear();
  zlhist_metJetUp.clear();
  zlhist_metJetDw.clear();
  zlhist_metResUp.clear();
  zlhist_metResDw.clear();
  zlhist_metUncUp.clear();
  zlhist_metUncDw.clear();
  tthist.clear();
  tthist_metJetUp.clear();
  tthist_metJetDw.clear();
  tthist_metResUp.clear();
  tthist_metResDw.clear();
  tthist_metUncUp.clear();
  tthist_metUncDw.clear();
  tthist_alt.clear();
  tthist_alt_metJetUp.clear();
  tthist_alt_metJetDw.clear();
  tthist_alt_metResUp.clear();
  tthist_alt_metResDw.clear();
  tthist_alt_metUncUp.clear();
  tthist_alt_metUncDw.clear();
  dihist.clear();
  dihist_metJetUp.clear();
  dihist_metJetDw.clear();
  dihist_metResUp.clear();
  dihist_metResDw.clear();
  dihist_metUncUp.clear();
  dihist_metUncDw.clear();
  gmhist.clear();
  gmhist_metJetUp.clear();
  gmhist_metJetDw.clear();
  gmhist_metResUp.clear();
  gmhist_metResDw.clear();
  gmhist_metUncUp.clear();
  gmhist_metUncDw.clear();
  qcdhist.clear();
  ewkwhist.clear();
  ewkwlhist_mu.clear();
  ewkwlhist_el.clear();
  ewkwlhist_ta.clear();
  ewkzhist.clear();
  vghist.clear();
  dthist.clear();

  for(auto hist_2D: znhist_2D) delete hist_2D;
  for(auto hist_2D: wlhist_2D) delete hist_2D;
  for(auto hist_2D: zlhist_2D) delete hist_2D;
  for(auto hist_2D: zlhist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D: zlhist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D: zlhist_metResUp_2D) delete hist_2D;
  for(auto hist_2D: zlhist_metResDw_2D) delete hist_2D;
  for(auto hist_2D: zlhist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D: zlhist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D: tthist_2D) delete hist_2D;
  for(auto hist_2D: tthist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D: tthist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D: tthist_metResUp_2D) delete hist_2D;
  for(auto hist_2D: tthist_metResDw_2D) delete hist_2D;
  for(auto hist_2D: tthist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D: tthist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_metJetUp_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_metJetDw_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_metResUp_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_metResDw_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_metUncUp_2D) delete hist_2D;
  for(auto hist_2D: tthist_alt_metUncDw_2D) delete hist_2D;
  for(auto hist_2D: dihist_2D) delete hist_2D;
  for(auto hist_2D: dihist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D: dihist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D: dihist_metResUp_2D) delete hist_2D;
  for(auto hist_2D: dihist_metResDw_2D) delete hist_2D;
  for(auto hist_2D: dihist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D: dihist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D: gmhist_2D) delete hist_2D;
  for(auto hist_2D: gmhist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D: gmhist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D: gmhist_metResUp_2D) delete hist_2D;
  for(auto hist_2D: gmhist_metResDw_2D) delete hist_2D;
  for(auto hist_2D: gmhist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D: gmhist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D: qcdhist_2D) delete hist_2D;
  for(auto hist_2D: ewkwhist_2D) delete hist_2D;
  for(auto hist_2D: ewkzhist_2D) delete hist_2D;
  for(auto hist_2D: vghist_2D) delete hist_2D;
  for(auto hist_2D: dthist_2D) delete hist_2D;

  znhist_2D.clear();
  wlhist_2D.clear();
  wlhist_mu_2D.clear();
  wlhist_el_2D.clear();
  wlhist_ta_2D.clear();
  zlhist_2D.clear();
  zlhist_metJetUp_2D.clear();
  zlhist_metJetDw_2D.clear();
  zlhist_metResUp_2D.clear();
  zlhist_metResDw_2D.clear();
  zlhist_metUncUp_2D.clear();
  zlhist_metUncDw_2D.clear();
  tthist_2D.clear();
  tthist_metJetUp_2D.clear();
  tthist_metJetDw_2D.clear();
  tthist_metResUp_2D.clear();
  tthist_metResDw_2D.clear();
  tthist_metUncUp_2D.clear();
  tthist_metUncDw_2D.clear();
  tthist_alt_2D.clear();
  tthist_alt_metJetUp_2D.clear();
  tthist_alt_metJetDw_2D.clear();
  tthist_alt_metResUp_2D.clear();
  tthist_alt_metResDw_2D.clear();
  tthist_alt_metUncUp_2D.clear();
  tthist_alt_metUncDw_2D.clear();
  dihist_2D.clear();
  dihist_metJetUp_2D.clear();
  dihist_metJetDw_2D.clear();
  dihist_metResUp_2D.clear();
  dihist_metResDw_2D.clear();
  dihist_metUncUp_2D.clear();
  dihist_metUncDw_2D.clear();
  gmhist_2D.clear();
  gmhist_metJetUp_2D.clear();
  gmhist_metJetDw_2D.clear();
  gmhist_metResUp_2D.clear();
  gmhist_metResDw_2D.clear();
  gmhist_metUncUp_2D.clear();
  gmhist_metUncDw_2D.clear();
  qcdhist_2D.clear();
  ewkwhist_2D.clear();
  ewkwlhist_mu_2D.clear();
  ewkwlhist_el_2D.clear();
  ewkwlhist_ta_2D.clear();
  ewkzhist_2D.clear();
  vghist_2D.clear();
  dthist_2D.clear();

  if(kfactzjet_vbf)
    kfactzjet_vbf->Close();
  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();
  if(kfactgjet_vbf)
    kfactgjet_vbf->Close();
  if(kffile)
    kffile->Close();
  if(kffile_alt)
    kffile_alt->Close();
  if(kffile_alt2)
    kffile_alt2->Close();
  if(kffile_unlops)
    kffile_unlops->Close();
  if(kffile_zewk)
    kffile_zewk->Close();
  if(kffile_wewk)
    kffile_wewk->Close();

  if(zntree) delete zntree;
  if(wltree) delete wltree;
  if(zltree) delete zltree;
  if(tttree) delete tttree;
  if(tttree_alt) delete tttree_alt;
  if(ditree) delete ditree;
  if(gmtree) delete gmtree;
  if(qcdtree) delete qcdtree;
  if(ewkwtree) delete ewkwtree;
  if(ewkztree) delete ewkztree;
  if(vgamtree) delete vgamtree;
  if(zgtree) delete zgtree;
  if(wgtree) delete wgtree;
  if(dttree) delete dttree;

  cout << "Templates for the signal region computed ..." << endl;
}

// build templates for photon+jets control region                                                                                                                           
void gamdatamchist(TFile* outfile,
		   const Category & category,
                   vector<string> observables,
                   vector<string> observables_2D,
		   const SamplesNLO & nloSamples,
                   const double & lumi   = 12.9,
		   const bool & isHInv   = false,
		   const bool & useJetHT = false,
		   const bool & applyPFWeight = false,
		   const bool & useTheoriestKFactors = false		   
                   ) {


  TChain* dttree   = new TChain("tree/tree");
  TChain* gmtree = new TChain("tree/tree");
  TChain* wgtree = new TChain("tree/tree");
  TChain* zgtree = new TChain("tree/tree");
  TChain* vltree = new TChain("tree/tree");

  //dttree->Add((baseInputTreePath+"/SinglePhoton_jecReReco/gamfilter/*root").c_str());
  dttree->Add((baseInputTreePath+"/SinglePhoton_reMiniAOD/gamfilter/*root").c_str());
  if(useJetHT)
    dttree->Add((baseInputTreePath+"/JetHT_jecReReco/gamfilter/*root").c_str());
  
  gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/*root").c_str());
  zgtree->Add((baseInputTreePath+"/ZnunuGJets_highPt/gamfilter/*root").c_str());
  wgtree->Add((baseInputTreePath+"/WGJets_highPt/gamfilter/*root").c_str());
  vltree->Add((baseInputTreePath+"/WJets/gamfilter/*root").c_str());

  vector<TH1*> dthist;
  vector<TH1*> qcdhist;
  vector<TH1*> gmhist;
  vector<TH1*> vghist;
  vector<TH1*> vlhist;

  vector<TH2*> dthist_2D;
  vector<TH2*> qcdhist_2D;
  vector<TH2*> vghist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> gmhist_2D;

  vector<double> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;

    TH1F* gmhist_temp = new TH1F(("gbkghistgam_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* qchist_temp = new TH1F(("qbkghistgam_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vghist_temp = new TH1F(("vgbkghistgam_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vlhist_temp = new TH1F(("vlbkghistgam_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dthist_temp = new TH1F(("datahistgam_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

    qcdhist.push_back(dynamic_cast<TH1*>(qchist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    vghist.push_back(dynamic_cast<TH1*>(vghist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
  }

 
  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty() )
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;

    TH2F* gmhist_temp = new TH2F(("gbkghistgam_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* qchist_temp = new TH2F(("qbkghistgam_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vghist_temp = new TH2F(("vgbkghistgam_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vlhist_temp = new TH2F(("vlbkghistgam_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dthist_temp = new TH2F(("datahistgam_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

    qcdhist_2D.push_back(dynamic_cast<TH2*>(qchist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vghist_2D.push_back(dynamic_cast<TH2*>(vghist_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));
  }


  // k-factors file from generator lebel: Z-boson pt at LO,NLO QCD and NLO QCD+EWK                                                                                           
  // get k-factors NLO                                                                                                                                                        

  TFile* kffile = NULL;
  TFile* kffile_alt = NULL;
  TFile* kffile_unlops = NULL;
  TH1*  alohist  = NULL;
  TH1*  anlohist = NULL;
  TH1*  aunlopshist = NULL;
  TH1*  aewkhist = NULL;
  TH1*  wnlohist = NULL;
  TH1*  wlohist  = NULL;
  TH1*  wewkhist  = NULL;
  TH1*  reweight_wln = NULL;

  if(not useTheoriestKFactors){
    kffile    = TFile::Open(kfactorFile.c_str());
    kffile_unlops = TFile::Open(kfactorFileUNLOPS.c_str());
    aunlopshist = (TH1*) kffile_unlops->Get("Func");
    alohist   = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    anlohist  = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    aewkhist  = (TH1*) kffile->Get("EWKcorr/photon");
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);

    // take open loop hitogram                                                                                                                                                    
    wnlohist = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
    wlohist  = (TH1*) kffile->Get("WJets_LO/inv_pt");
    wewkhist = (TH1*) kffile->Get("EWKcorr/W");

    if(wewkhist)
      wewkhist->Divide(wnlohist);
    if(wnlohist)
      wnlohist->Divide(wlohist);
  }
  else{

    kffile_alt = TFile::Open(kFactorTheoristFile_wln.c_str());
    TH1* kfact_nloqcd_wln = (TH1*) kffile_alt->Get("evj_pTV_K_NLO");
    TH1* kfact_nloewk_wln = (TH1*) kffile_alt->Get("evj_pTV_kappa_NLO_EW");
    TH1* kfact_sudewk_wln = (TH1*) kffile_alt->Get("evj_pTV_kappa_NNLO_Sud");
    reweight_wln = (TH1*) kfact_nloqcd_wln->Clone("reweight_wln");
    reweight_wln->Reset("ICES");
    
    for(int iBin = 1; iBin < reweight_wln->GetNbinsX(); iBin++)
      reweight_wln->SetBinContent(iBin,kfact_nloqcd_wln->GetBinContent(iBin)*(1+kfact_nloewk_wln->GetBinContent(iBin)+kfact_sudewk_wln->GetBinContent(iBin)));
    
    // no gamma+jets k-factors yet
    anlohist  = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    alohist  = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    aewkhist  = (TH1*) kffile->Get("EWKcorr/photon");
    
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);

  }

  vector<TH1*> ahists;
  vector<TH1*> whists;
  vector<TH1*> ehists;

  if(not useTheoriestKFactors){
    if(nloSamples.usePhotonJetsNLO){
      ahists.push_back(aewkhist);
      ahists.push_back(aunlopshist);
    }
    else{
      ahists.push_back(aewkhist);
      ahists.push_back(anlohist);
      ahists.push_back(aunlopshist);
    }
    if(nloSamples.useWJetsNLO){
      whists.push_back(wewkhist);
    }
    else{
      whists.push_back(wnlohist);
      whists.push_back(wewkhist);
    }
  }
  else{
    ahists.push_back(aewkhist);
    ahists.push_back(anlohist);    
    whists.push_back(reweight_wln);
  }

  TFile* kfactwjet_vbf = NULL;
  TFile* kfactgjet_vbf = NULL;

  if((category == Category::VBF or category == Category::VBFrelaxed)  and not useTheoriestKFactors){ // apply further k-factors going to the VBF selections                         
                                                                     
    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets_v2.root");
    kfactgjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_gjets_v2.root");

    TH1F* wjet_nlo_vbf = (TH1F*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      wjet_nlo_vbf = (TH1F*) kfactwjet_vbf->Get("bosonPt_NLO_vbf_relaxed");
    TH1F* wjet_nlo_mj  = (TH1F*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");

    if(category == Category::VBF)
      wjet_nlo_vbf->Divide((TH1F*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      wjet_nlo_vbf->Divide((TH1F*) kfactwjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
    
    wjet_nlo_mj->Divide((TH1F*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);

    if(not nloSamples.useWJetsNLO)
      whists.push_back(wjet_nlo_vbf);

    TH1F* gjet_nlo_vbf = (TH1F*) kfactgjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      gjet_nlo_vbf = (TH1F*) kfactgjet_vbf->Get("bosonPt_NLO_vbf_relaxed");

    TH1F* gjet_nlo_mj  = (TH1F*) kfactgjet_vbf->Get("bosonPt_NLO_monojet");
    if(category == Category::VBF)
      gjet_nlo_vbf->Divide((TH1F*) kfactgjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      gjet_nlo_vbf->Divide((TH1F*) kfactgjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
    
    gjet_nlo_mj->Divide((TH1F*) kfactgjet_vbf->Get("bosonPt_LO_monojet"));
    gjet_nlo_vbf->Divide(gjet_nlo_mj);
    if(not nloSamples.usePhotonJetsNLO)
      ahists.push_back(gjet_nlo_vbf);
  }
  
  // NLO k-factor for Zgamma is a function of gamma pT
  cout<<"gamma+jets control region --> gamma+jets"<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,Sample::gam,category,false,1.00,lumi,ahists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"gamma+jets control region --> Wgamma+jets"<<endl;
  // NLO k-factor for Wgamma is a flat 34% from photon pT > 175 GeV
  makehist4(wgtree,vghist,vghist_2D,true,Sample::gam,category,false,1.34,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"gamma+jets control region --> Zgamma+jets"<<endl;
  makehist4(zgtree,vghist,vghist_2D,true,Sample::gam,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"gamma+jets control region --> QCD"<<endl;
  makehist4(dttree,qcdhist,qcdhist_2D,false,Sample::qcdgam,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"gamma+jets control region --> W+jets"<<endl;
  makehist4(vltree,vlhist,vlhist_2D,true,Sample::gam,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"gamma+jets control region --> data"<<endl;
  makehist4(dttree,dthist,dthist_2D,false,Sample::gam,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  outfile->cd();
  if(not outfile->GetDirectory("GJ"))
    outfile->mkdir("GJ");
  outfile->cd("GJ");

  for(auto hist : dthist) hist->Write();
  for(auto hist : gmhist) hist->Write();
  for(auto hist : vghist) hist->Write();
  for(auto hist : vlhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();

  for(auto hist : dthist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
  for(auto hist : gmhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
  for(auto hist : vghist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
  for(auto hist : vlhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
  for(auto hist : qcdhist_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

  outfile->cd();

  whists.clear();
  ahists.clear();

  for(auto hist : dthist) delete hist;
  for(auto hist : qcdhist) delete hist;
  for(auto hist : vghist) delete hist;
  for(auto hist : vlhist) delete hist;
  for(auto hist : gmhist) delete hist;

  for(auto hist : dthist_2D) delete hist;
  for(auto hist : qcdhist_2D) delete hist;
  for(auto hist : vghist_2D) delete hist;
  for(auto hist : vlhist_2D) delete hist;
  for(auto hist : gmhist_2D) delete hist;

  dthist.clear();
  qcdhist.clear();
  vghist.clear();
  vlhist.clear();
  gmhist.clear();
  dthist_2D.clear();
  qcdhist_2D.clear();
  vghist_2D.clear();
  vlhist_2D.clear();
  gmhist_2D.clear();

  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();
  if(kfactgjet_vbf)
    kfactgjet_vbf->Close();
  if(kffile)
    kffile->Close();
  if(kffile_alt)
    kffile->Close();
  if(kffile_unlops)
    kffile_unlops->Close();
  
  if(dttree) delete dttree;
  if(gmtree) delete gmtree;
  if(wgtree) delete wgtree;
  if(zgtree) delete zgtree;
  if(vltree) delete vltree;

  cout << "Templates for the gamma+jets control region computed ..." << endl;
}

//build templates for Zmumu,Zee,Wenu,Wmunu                                                                                                                                  
void lepdatamchist(TFile* outfile,
		   const Sample & sample,
		   const Category & category,
		   vector<string> observables,
		   vector<string> observables_2D,
		   const double & lumi,
		   const SamplesNLO & nloSamples,
		   const bool & doShapeSystematics = false,
		   const bool & isHInv             = false,
		   const bool & useSinglePhoton    = false,
		   const bool & useJetHT           = false,
		   const bool & applyPFWeight      = false,
		   const bool & useTheoriestKFactors = false) {

  if (sample != Sample::zmm && sample != Sample::zee && sample != Sample::wmn && sample != Sample::wen) return;

  if(useJetHT and useSinglePhoton){
    cerr<<"Error: useJetHT and useSinglePhoton can't be both true -- return "<<endl;
    return;
  }

  TChain* tttree  = new TChain("tree/tree");
  TChain* dbtree  = new TChain("tree/tree"); 
  TChain* gmtree  = new TChain("tree/tree");
  TChain* qctree  = new TChain("tree/tree");
  TChain* vltree  = new TChain("tree/tree");
  TChain* vlltree = new TChain("tree/tree");
  TChain* ewkwtree = new TChain("tree/tree");
  TChain* ewkztree = new TChain("tree/tree");
  TChain* zgtree  = new TChain("tree/tree");
  TChain* wgtree  = new TChain("tree/tree");

  TChain* vltree_nlo1  = new TChain("tree/tree");
  TChain* vltree_nlo2  = new TChain("tree/tree");
  TChain* vltree_nlo3  = new TChain("tree/tree");
  TChain* vltree_nlo4  = new TChain("tree/tree");

  TChain* vlltree_nlo1  = new TChain("tree/tree");
  TChain* vlltree_nlo2  = new TChain("tree/tree");
  TChain* vlltree_nlo3  = new TChain("tree/tree");
  TChain* vlltree_nlo4  = new TChain("tree/tree");

  TChain* dttree  = new TChain("tree/tree");
  TChain* dttree_2  = NULL;

  string suffix;
  if(sample == Sample::zmm){

    suffix = "zmm";
    if(nloSamples.useDYJetsNLO){
      vlltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/*Pt-100To250*root").c_str());
      vlltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/*Pt-250To400*root").c_str());
      vlltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/*Pt-400To650*root").c_str());
      vlltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/*Pt-650ToInf*root").c_str());
    }
    else{
      vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/*root").c_str());
    }


    if(nloSamples.useWJetsNLO){
      vltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zmmfilter/*Pt-100To250*root").c_str());
      vltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zmmfilter/*Pt-250To400*root").c_str());
      vltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zmmfilter/*Pt-400To600*root").c_str());
      vltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zmmfilter/*Pt-600ToInf*root").c_str());
    }
    else{
      vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zmmfilter/*root").c_str());
    }

    gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/zmmfilter/*root").c_str());
    qctree->Add((baseInputTreePath+"/QCD/zmmfilter/*root").c_str());
    dbtree->Add((baseInputTreePath+"/DiBoson/zmmfilter/*root").c_str());
    tttree->Add((baseInputTreePath+"/TopAlternative/zmmfilter/*root").c_str());
    ewkwtree->Add((baseInputTreePath+"/WJetsEWK/zmmfilter/*root").c_str());
    ewkztree->Add((baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/*root").c_str());
    zgtree->Add((baseInputTreePath+"/ZLLGJets/zmmfilter/*root").c_str());
    wgtree->Add((baseInputTreePath+"/WGJets/zmmfilter/*root").c_str());
    dttree->Add((baseInputTreePath+"/MET_reMiniAOD/zmmfilter/*root").c_str());
  }
  else if(sample == Sample::wmn){

    suffix = "wmn";
    if(nloSamples.useDYJetsNLO){
      vlltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wmnfilter/*Pt-100To250*root").c_str());
      vlltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wmnfilter/*Pt-250To400*root").c_str());
      vlltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wmnfilter/*Pt-400To650*root").c_str());
      vlltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wmnfilter/*Pt-650ToInf*root").c_str());
    }
    else{
      vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wmnfilter/*root").c_str());
    }


    if(nloSamples.useWJetsNLO){
      vltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/*Pt-100To250*root").c_str());
      vltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/*Pt-250To400*root").c_str());
      vltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/*Pt-400To600*root").c_str());
      vltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/*Pt-600ToInf*root").c_str());
    }
    else{
      vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/*root").c_str());
    }

    qctree->Add((baseInputTreePath+"/QCD/wmnfilter/*root").c_str());
    dbtree->Add((baseInputTreePath+"/DiBoson/wmnfilter/*root").c_str());
    gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/wmnfilter/*root").c_str());
    tttree->Add((baseInputTreePath+"/TopAlternative/wmnfilter/*root").c_str());
    ewkwtree->Add((baseInputTreePath+"/WJetsEWK/wmnfilter/*root").c_str());
    ewkztree->Add((baseInputTreePath+"/ZJetsToLLEWK/wmnfilter/*root").c_str());
    zgtree->Add((baseInputTreePath+"/ZLLGJets/wmnfilter/*root").c_str());
    wgtree->Add((baseInputTreePath+"/WGJets/wmnfilter/*root").c_str());
    dttree->Add((baseInputTreePath+"/MET_reMiniAOD/wmnfilter/*root").c_str());
  }
  else if(sample == Sample::zee){

    suffix = "zee";
    if(nloSamples.useDYJetsNLO){
      vlltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/*Pt-100To250*root").c_str());
      vlltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/*Pt-250To400*root").c_str());
      vlltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/*Pt-400To650*root").c_str());
      vlltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/*Pt-650ToInf*root").c_str());
    }
    else{
      vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/*root").c_str());
    }

    if(nloSamples.useWJetsNLO){
      vltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zeefilter/*Pt-100To250*root").c_str());
      vltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zeefilter/*Pt-250To400*root").c_str());
      vltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zeefilter/*Pt-400To600*root").c_str());
      vltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zeefilter/*Pt-600ToInf*root").c_str());
    }
    else{
      vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/zeefilter/*root").c_str());
    }

    qctree->Add((baseInputTreePath+"/QCD/zeefilter/*root").c_str());
    dbtree->Add((baseInputTreePath+"/DiBoson/zeefilter/*root").c_str());
    gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/zeefilter/*root").c_str());
    tttree->Add((baseInputTreePath+"/TopAlternative/zeefilter/*root").c_str());
    ewkwtree->Add((baseInputTreePath+"/WJetsEWK/zeefilter/*root").c_str());
    ewkztree->Add((baseInputTreePath+"/ZJetsToLLEWK/zeefilter/*root").c_str());
    zgtree->Add((baseInputTreePath+"/ZLLGJets/zeefilter/*root").c_str());
    wgtree->Add((baseInputTreePath+"/WGJets/zeefilter/*root").c_str());
    dttree->Add((baseInputTreePath+"/SingleElectron_reMiniAOD/zeefilter/*root").c_str());
    dttree_2 = new TChain("tree/tree");
    if(useSinglePhoton)
      dttree_2->Add((baseInputTreePath+"/SinglePhoton_jecReReco/zeefilter/*root").c_str());
    else if(useJetHT)
      dttree_2->Add((baseInputTreePath+"/JetHT_jecReReco/zeefilter/*root").c_str());
  }
  else if(sample == Sample::wen){

    suffix = "wen";

    if(nloSamples.useDYJetsNLO){
      vlltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wenfilter/*Pt-100To250*root").c_str());
      vlltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wenfilter/*Pt-250To400*root").c_str());
      vlltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wenfilter/*Pt-400To650*root").c_str());
      vlltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wenfilter/*Pt-650ToInf*root").c_str());
    }
    else{
      vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/wenfilter/*root").c_str());
    }


    if(nloSamples.useWJetsNLO){
      vltree_nlo1->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/*Pt-100To250*root").c_str());
      vltree_nlo2->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/*Pt-250To400*root").c_str());
      vltree_nlo3->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/*Pt-400To600*root").c_str());
      vltree_nlo4->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/*Pt-600ToInf*root").c_str());
    }
    else{
      vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/*root").c_str());
    }

    qctree->Add((baseInputTreePath+"/QCD/wenfilter/*root").c_str());
    dbtree->Add((baseInputTreePath+"/DiBoson/wenfilter/*root").c_str());
    gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/wenfilter/*root").c_str());
    tttree->Add((baseInputTreePath+"/TopAlternative/wenfilter/*root").c_str());

    ewkwtree->Add((baseInputTreePath+"/WJetsEWK/wenfilter/*root").c_str());
    ewkztree->Add((baseInputTreePath+"/ZJetsToLLEWK/wenfilter/*root").c_str());
    zgtree->Add((baseInputTreePath+"/ZLLGJets/wenfilter/*root").c_str());
    wgtree->Add((baseInputTreePath+"/WGJets/wenfilter/*root").c_str());

    dttree->Add((baseInputTreePath+"/SingleElectron_reMiniAOD/wenfilter/*root").c_str());
    dttree_2 = new TChain("tree/tree");
    if(useJetHT)
      dttree_2->Add((baseInputTreePath+"/JetHT_jecReReco/wenfilter/*root").c_str());
    else if(useSinglePhoton)
      dttree_2->Add((baseInputTreePath+"/SinglePhoton_jecReReco/wenfilter/*root").c_str());
  }

  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_metJetUp;
  vector<TH1*> tthist_metJetDw;
  vector<TH1*> tthist_metResUp;
  vector<TH1*> tthist_metResDw;
  vector<TH1*> tthist_metUncUp;
  vector<TH1*> tthist_metUncDw;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> dbhist_metJetUp;
  vector<TH1*> dbhist_metJetDw;
  vector<TH1*> dbhist_metResUp;
  vector<TH1*> dbhist_metResDw;
  vector<TH1*> dbhist_metUncUp;
  vector<TH1*> dbhist_metUncDw;
  vector<TH1*> gmhist;
  vector<TH1*> gmhist_metJetUp;
  vector<TH1*> gmhist_metJetDw;
  vector<TH1*> gmhist_metResUp;
  vector<TH1*> gmhist_metResDw;
  vector<TH1*> gmhist_metUncUp;
  vector<TH1*> gmhist_metUncDw;
  vector<TH1*> vlhist;
  vector<TH1*> vlhist_metJetUp;
  vector<TH1*> vlhist_metJetDw;
  vector<TH1*> vlhist_metResUp;
  vector<TH1*> vlhist_metResDw;
  vector<TH1*> vlhist_metUncUp;
  vector<TH1*> vlhist_metUncDw;
  vector<TH1*> vllhist;
  vector<TH1*> vllhist_metJetUp;
  vector<TH1*> vllhist_metJetDw;
  vector<TH1*> vllhist_metResUp;
  vector<TH1*> vllhist_metResDw;
  vector<TH1*> vllhist_metUncUp;
  vector<TH1*> vllhist_metUncDw;

  vector<TH1*> ewkwhist;
  vector<TH1*> ewkzhist;
  vector<TH1*> vghist;

  vector<double> bins;
  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;


    TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* gmhist_temp = new TH1F((string("gbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkwhist_temp = new TH1F((string("ewkwbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkzhist_temp = new TH1F((string("ewkzbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vghist_temp = new TH1F((string("vgbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));
    ewkwhist.push_back(dynamic_cast<TH1*>(ewkwhist_temp));
    ewkzhist.push_back(dynamic_cast<TH1*>(ewkzhist_temp));
    vghist.push_back(dynamic_cast<TH1*>(vghist_temp));

    if(doShapeSystematics){

      TH1F* tthist_metJetUp_temp = new TH1F((string("tbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metJetDw_temp = new TH1F((string("tbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metResUp_temp = new TH1F((string("tbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metResDw_temp = new TH1F((string("tbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metUncUp_temp = new TH1F((string("tbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_metUncDw_temp = new TH1F((string("tbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      tthist_metJetUp.push_back(dynamic_cast<TH1*>(tthist_metJetUp_temp));
      tthist_metJetDw.push_back(dynamic_cast<TH1*>(tthist_metJetDw_temp));
      tthist_metResUp.push_back(dynamic_cast<TH1*>(tthist_metResUp_temp));
      tthist_metResDw.push_back(dynamic_cast<TH1*>(tthist_metResDw_temp));
      tthist_metUncUp.push_back(dynamic_cast<TH1*>(tthist_metUncUp_temp));
      tthist_metUncDw.push_back(dynamic_cast<TH1*>(tthist_metUncDw_temp));

      TH1F* dbhist_metJetUp_temp = new TH1F((string("dbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metJetDw_temp = new TH1F((string("dbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metResUp_temp = new TH1F((string("dbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metResDw_temp = new TH1F((string("dbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metUncUp_temp = new TH1F((string("dbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metUncDw_temp = new TH1F((string("dbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      dbhist_metJetUp.push_back(dynamic_cast<TH1*>(dbhist_metJetUp_temp));
      dbhist_metJetDw.push_back(dynamic_cast<TH1*>(dbhist_metJetDw_temp));
      dbhist_metResUp.push_back(dynamic_cast<TH1*>(dbhist_metResUp_temp));
      dbhist_metResDw.push_back(dynamic_cast<TH1*>(dbhist_metResDw_temp));
      dbhist_metUncUp.push_back(dynamic_cast<TH1*>(dbhist_metUncUp_temp));
      dbhist_metUncDw.push_back(dynamic_cast<TH1*>(dbhist_metUncDw_temp));

      TH1F* gmhist_metJetUp_temp = new TH1F((string("gbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metJetDw_temp = new TH1F((string("gbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metResUp_temp = new TH1F((string("gbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metResDw_temp = new TH1F((string("gbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metUncUp_temp = new TH1F((string("gbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metUncDw_temp = new TH1F((string("gbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      gmhist_metJetUp.push_back(dynamic_cast<TH1*>(gmhist_metJetUp_temp));
      gmhist_metJetDw.push_back(dynamic_cast<TH1*>(gmhist_metJetDw_temp));
      gmhist_metResUp.push_back(dynamic_cast<TH1*>(gmhist_metResUp_temp));
      gmhist_metResDw.push_back(dynamic_cast<TH1*>(gmhist_metResDw_temp));
      gmhist_metUncUp.push_back(dynamic_cast<TH1*>(gmhist_metUncUp_temp));
      gmhist_metUncDw.push_back(dynamic_cast<TH1*>(gmhist_metUncDw_temp));

    }

    if((sample == Sample::zmm or sample == Sample::zee) and doShapeSystematics){

      TH1F* vlhist_metJetUp_temp = new TH1F((string("vlbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metJetDw_temp = new TH1F((string("vlbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metResUp_temp = new TH1F((string("vlbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metResDw_temp = new TH1F((string("vlbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metUncUp_temp = new TH1F((string("vlbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metUncDw_temp = new TH1F((string("vlbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      vlhist_metJetUp.push_back(dynamic_cast<TH1*>(vlhist_metJetUp_temp));
      vlhist_metJetDw.push_back(dynamic_cast<TH1*>(vlhist_metJetDw_temp));
      vlhist_metResUp.push_back(dynamic_cast<TH1*>(vlhist_metResUp_temp));
      vlhist_metResDw.push_back(dynamic_cast<TH1*>(vlhist_metResDw_temp));
      vlhist_metUncUp.push_back(dynamic_cast<TH1*>(vlhist_metUncUp_temp));
      vlhist_metUncDw.push_back(dynamic_cast<TH1*>(vlhist_metUncDw_temp));
      
    }
    else if((sample == Sample::wmn or sample == Sample::wen) and doShapeSystematics){

      TH1F* vllhist_metJetUp_temp = new TH1F((string("vllbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metJetDw_temp = new TH1F((string("vllbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metResUp_temp = new TH1F((string("vllbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metResDw_temp = new TH1F((string("vllbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metUncUp_temp = new TH1F((string("vllbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metUncDw_temp = new TH1F((string("vllbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

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
  vector<TH2*> tthist_metJetUp_2D;
  vector<TH2*> tthist_metJetDw_2D;
  vector<TH2*> tthist_metResUp_2D;
  vector<TH2*> tthist_metResDw_2D;
  vector<TH2*> tthist_metUncUp_2D;
  vector<TH2*> tthist_metUncDw_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> dbhist_metJetUp_2D;
  vector<TH2*> dbhist_metJetDw_2D;
  vector<TH2*> dbhist_metResUp_2D;
  vector<TH2*> dbhist_metResDw_2D;
  vector<TH2*> dbhist_metUncUp_2D;
  vector<TH2*> dbhist_metUncDw_2D; 
  vector<TH2*> gmhist_2D;
  vector<TH2*> gmhist_metJetUp_2D;
  vector<TH2*> gmhist_metJetDw_2D;
  vector<TH2*> gmhist_metResUp_2D;
  vector<TH2*> gmhist_metResDw_2D;
  vector<TH2*> gmhist_metUncUp_2D;
  vector<TH2*> gmhist_metUncDw_2D; 
  vector<TH2*> vlhist_2D; 
  vector<TH2*> vlhist_metJetUp_2D;
  vector<TH2*> vlhist_metJetDw_2D;
  vector<TH2*> vlhist_metResUp_2D;
  vector<TH2*> vlhist_metResDw_2D;
  vector<TH2*> vlhist_metUncUp_2D;
  vector<TH2*> vlhist_metUncDw_2D;
  vector<TH2*> vllhist_2D;
  vector<TH2*> vllhist_metJetUp_2D;
  vector<TH2*> vllhist_metJetDw_2D;
  vector<TH2*> vllhist_metResUp_2D;
  vector<TH2*> vllhist_metResDw_2D;
  vector<TH2*> vllhist_metUncUp_2D;
  vector<TH2*> vllhist_metUncDw_2D;

  vector<TH2*> ewkwhist_2D;
  vector<TH2*> ewkzhist_2D;
  vector<TH2*> vghist_2D;

  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;


    TH2F* dthist_temp  = new TH2F((string("datahist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* tthist_temp  = new TH2F((string("tbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dbhist_temp  = new TH2F((string("dbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* gmhist_temp  = new TH2F((string("gbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* qchist_temp  = new TH2F((string("qbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vlhist_temp  = new TH2F((string("vlbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vllhist_temp = new TH2F((string("vllbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkwhist_temp  = new TH2F((string("ewkwbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkzhist_temp = new TH2F((string("ewkzbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vghist_temp = new TH2F((string("vgbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    dbhist_2D.push_back(dynamic_cast<TH2*>(dbhist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    qchist_2D.push_back(dynamic_cast<TH2*>(qchist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vllhist_2D.push_back(dynamic_cast<TH2*>(vllhist_temp));
    ewkwhist_2D.push_back(dynamic_cast<TH2*>(ewkwhist_temp));
    ewkzhist_2D.push_back(dynamic_cast<TH2*>(ewkzhist_temp));
    vghist_2D.push_back(dynamic_cast<TH2*>(vghist_temp));
    
    if(doShapeSystematics){

      TH2F* tthist_metJetUp_temp = new TH2F((string("tbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metJetDw_temp = new TH2F((string("tbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metResUp_temp = new TH2F((string("tbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metResDw_temp = new TH2F((string("tbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metUncUp_temp = new TH2F((string("tbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* tthist_metUncDw_temp = new TH2F((string("tbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

      tthist_metJetUp_2D.push_back(dynamic_cast<TH2*>(tthist_metJetUp_temp));
      tthist_metJetDw_2D.push_back(dynamic_cast<TH2*>(tthist_metJetDw_temp));
      tthist_metResUp_2D.push_back(dynamic_cast<TH2*>(tthist_metResUp_temp));
      tthist_metResDw_2D.push_back(dynamic_cast<TH2*>(tthist_metResDw_temp));
      tthist_metUncUp_2D.push_back(dynamic_cast<TH2*>(tthist_metUncUp_temp));
      tthist_metUncDw_2D.push_back(dynamic_cast<TH2*>(tthist_metUncDw_temp));

      TH2F* dbhist_metJetUp_temp = new TH2F((string("dbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metJetDw_temp = new TH2F((string("dbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metResUp_temp = new TH2F((string("dbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metResDw_temp = new TH2F((string("dbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metUncUp_temp = new TH2F((string("dbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metUncDw_temp = new TH2F((string("dbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

      dbhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetUp_temp));
      dbhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetDw_temp));
      dbhist_metResUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metResUp_temp));
      dbhist_metResDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metResDw_temp));
      dbhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncUp_temp));
      dbhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncDw_temp));

      TH2F* gmhist_metJetUp_temp = new TH2F((string("gbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metJetDw_temp = new TH2F((string("gbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metResUp_temp = new TH2F((string("gbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metResDw_temp = new TH2F((string("gbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metUncUp_temp = new TH2F((string("gbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metUncDw_temp = new TH2F((string("gbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

      gmhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetUp_temp));
      gmhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetDw_temp));
      gmhist_metResUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metResUp_temp));
      gmhist_metResDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metResDw_temp));
      gmhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncUp_temp));
      gmhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncDw_temp));

    }

    if((sample == Sample::zmm or sample == Sample::zee) and doShapeSystematics){

      TH2F* vlhist_metJetUp_temp = new TH2F((string("vlbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metJetDw_temp = new TH2F((string("vlbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metResUp_temp = new TH2F((string("vlbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metResDw_temp = new TH2F((string("vlbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metUncUp_temp = new TH2F((string("vlbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metUncDw_temp = new TH2F((string("vlbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

      vlhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetUp_temp));
      vlhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetDw_temp));
      vlhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metResUp_temp));
      vlhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metResDw_temp));
      vlhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncUp_temp));
      vlhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncDw_temp));
      
    }
    else if((sample == Sample::wmn or sample == Sample::wen) and doShapeSystematics){

      TH2F* vllhist_metJetUp_temp = new TH2F((string("vllbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metJetDw_temp = new TH2F((string("vllbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metResUp_temp = new TH2F((string("vllbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metResDw_temp = new TH2F((string("vllbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metUncUp_temp = new TH2F((string("vllbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metUncDw_temp = new TH2F((string("vllbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

      vllhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetUp_temp));
      vllhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetDw_temp));
      vllhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metResUp_temp));
      vllhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metResDw_temp));
      vllhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncUp_temp));
      vllhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncDw_temp));
     
    }
  }

  // make a chain of two files
  TChain* dttree_merged = new TChain("tree");
  dttree_merged->Add(dttree);
  if(dttree_2)
    dttree_merged->Add(dttree_2);

  TFile* kffile = NULL;
  TFile* kffile_alt = NULL;
  TFile* kffile_unlops = NULL;
  TH1*  znlohist = NULL;
  TH1*  zlohist  = NULL;
  TH1* zewkhist  = NULL;
  TH1*  wnlohist = NULL;
  TH1*  wlohist  = NULL;
  TH1* wewkhist  = NULL;
  TH1* anlohist  = NULL;
  TH1*  alohist  = NULL;
  TH1* aewkhist  = NULL;
  TH1* aunlopshist = NULL;
  TH1* reweight_zvv =  NULL;
  TH1* reweight_zll =  NULL;
  TH1* reweight_wln =  NULL;

  if(not useTheoriestKFactors){

    kffile = TFile::Open(kfactorFile.c_str());
    znlohist = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
    zlohist  = (TH1*) kffile->Get("ZJets_LO/inv_pt");
    zewkhist = (TH1*) kffile->Get("EWKcorr/Z");

    if(zewkhist)
      zewkhist->Divide(znlohist);
    if(znlohist)
      znlohist->Divide(zlohist);

    wnlohist  = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
    wlohist   = (TH1*) kffile->Get("WJets_LO/inv_pt");
    wewkhist  = (TH1*) kffile->Get("EWKcorr/W");

    if(wewkhist)
      wewkhist->Divide(wnlohist);
    if(wnlohist)
      wnlohist->Divide(wlohist);
    
    kffile_unlops = TFile::Open(kfactorFileUNLOPS.c_str());
    anlohist  = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    alohist   = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    aewkhist  = (TH1*) kffile->Get("EWKcorr/photon");
    aunlopshist = (TH1*) kffile_unlops->Get("Func");
 
    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);
  }
  else{


    kffile = TFile::Open(kFactorTheoristFile_zll.c_str());
    TH1* kfact_nloqcd_zll = (TH1*) kffile->Get("eej_pTV_K_NLO");
    TH1* kfact_nloewk_zll = (TH1*) kffile->Get("eej_pTV_kappa_NLO_EW");
    TH1* kfact_sudewk_zll = (TH1*) kffile->Get("eej_pTV_kappa_NNLO_Sud");
    reweight_zll = (TH1*) kfact_nloqcd_zll->Clone("reweight_zll");
    reweight_zll->Reset("ICES");

    for(int iBin = 1; iBin < reweight_zll->GetNbinsX(); iBin++)
      reweight_zll->SetBinContent(iBin,kfact_nloqcd_zll->GetBinContent(iBin)*(1+kfact_nloewk_zll->GetBinContent(iBin)+kfact_sudewk_zll->GetBinContent(iBin)));

    kffile_alt = TFile::Open(kFactorTheoristFile_wln.c_str());
    TH1* kfact_nloqcd_wln = (TH1*) kffile_alt->Get("evj_pTV_K_NLO");
    TH1* kfact_nloewk_wln = (TH1*) kffile_alt->Get("evj_pTV_kappa_NLO_EW");
    TH1* kfact_sudewk_wln = (TH1*) kffile_alt->Get("evj_pTV_kappa_NNLO_Sud");
    reweight_wln = (TH1*) kfact_nloqcd_wln->Clone("reweight_wln");
    reweight_wln->Reset("ICES");

    for(int iBin = 1; iBin < reweight_wln->GetNbinsX(); iBin++)
      reweight_wln->SetBinContent(iBin,kfact_nloqcd_wln->GetBinContent(iBin)*(1+kfact_nloewk_wln->GetBinContent(iBin)+kfact_sudewk_wln->GetBinContent(iBin)));

    // no gamma+jets k-factors yet                                                                                                                                                                    
    anlohist  = (TH1*) kffile->Get("GJets_1j_NLO/nominal_G");
    alohist  = (TH1*) kffile->Get("GJets_LO/inv_pt_G");
    aewkhist  = (TH1*) kffile->Get("EWKcorr/photon");

    if(aewkhist)
      aewkhist->Divide(anlohist);
    if(anlohist)
      anlohist->Divide(alohist);
  }

  // take open loop hitogram                                                                                                                                                    
  vector<TH1*> ehists;
  vector<TH1*> vlhists;
  vector<TH1*> ahists;
  vector<TH1*> vllhists;

  if(not useTheoriestKFactors){
    // apply NLO QCD and EWK corrections for Zll and Wlnu                                                                                                                      
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);
    ahists.push_back(aunlopshist);
    vlhists.push_back(wnlohist);
    vlhists.push_back(wewkhist);
    vllhists.push_back(znlohist);
    vllhists.push_back(zewkhist);
    /////
    if(nloSamples.useWJetsNLO){ // no k-factors
      vlhists.clear();
      vlhists.push_back(wewkhist);
    }
    if(nloSamples.useDYJetsNLO){ // no k-factors
      vllhists.clear();
      vllhists.push_back(zewkhist);
    }
  }
  else{
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);    
    vlhists.push_back(reweight_wln);
    vllhists.push_back(reweight_zll);
  }

  // VBF k-factors
  TFile* kfactzjet_vbf = NULL;
  TFile* kfactwjet_vbf = NULL;
  TFile* kfactgjet_vbf = NULL;

  if((category == Category::VBF or category == Category::VBFrelaxed) and not useTheoriestKFactors){ // apply further k-factors going to the VBF selections                                                                                                                
    kfactzjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_zjets_v2.root");
    kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets_v2.root");
    kfactgjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_gjets_v2.root");

    TH1F* zjet_nlo_vbf = (TH1F*) kfactzjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      zjet_nlo_vbf = (TH1F*) kfactzjet_vbf->Get("bosonPt_NLO_vbf_relaxed");

    TH1F* zjet_nlo_mj  = (TH1F*) kfactzjet_vbf->Get("bosonPt_NLO_monojet");
    if(category == Category::VBF)
      zjet_nlo_vbf->Divide((TH1F*) kfactzjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      zjet_nlo_vbf->Divide((TH1F*) kfactzjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
    
    zjet_nlo_mj->Divide((TH1F*) kfactzjet_vbf->Get("bosonPt_LO_monojet"));
    zjet_nlo_vbf->Divide(zjet_nlo_mj);
    if(not nloSamples.useDYJetsNLO)
      vllhists.push_back(zjet_nlo_vbf);

    
    TH1F* wjet_nlo_vbf = (TH1F*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      wjet_nlo_vbf = (TH1F*) kfactwjet_vbf->Get("bosonPt_NLO_vbf_relaxed");

    TH1F* wjet_nlo_mj  = (TH1F*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
    if(category == Category::VBF)
      wjet_nlo_vbf->Divide((TH1F*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      wjet_nlo_vbf->Divide((TH1F*) kfactwjet_vbf->Get("bosonPt_LO_vbf_relaxed"));

    wjet_nlo_mj->Divide((TH1F*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));
    wjet_nlo_vbf->Divide(wjet_nlo_mj);
    if(not nloSamples.useWJetsNLO)
      vlhists.push_back(wjet_nlo_vbf);

    TH1F* gjet_nlo_vbf = (TH1F*) kfactgjet_vbf->Get("bosonPt_NLO_vbf");
    if(category == Category::VBFrelaxed)
      gjet_nlo_vbf = (TH1F*) kfactgjet_vbf->Get("bosonPt_NLO_vbf_relaxed");

    TH1F* gjet_nlo_mj  = (TH1F*) kfactgjet_vbf->Get("bosonPt_NLO_monojet");
    if(category == Category::VBF)
      gjet_nlo_vbf->Divide((TH1F*) kfactgjet_vbf->Get("bosonPt_LO_vbf"));
    else if(category == Category::VBFrelaxed)
      gjet_nlo_vbf->Divide((TH1F*) kfactgjet_vbf->Get("bosonPt_LO_vbf_relaxed"));
    
    gjet_nlo_mj->Divide((TH1F*) kfactgjet_vbf->Get("bosonPt_LO_monojet"));
    gjet_nlo_vbf->Divide(gjet_nlo_mj);
    if(not nloSamples.usePhotonJetsNLO)
      ahists.push_back(gjet_nlo_vbf);
  }
  
  /// v-EWK k-factord
  TFile* kffile_zewk = TFile::Open(kFactorFile_zjetewk.c_str(),"READ");
  TFile* kffile_wewk = TFile::Open(kFactorFile_wjetewk.c_str(),"READ");
  vector<TH2*> zewkhists;
  vector<TH2*> wewkhists;
  if(applyEWKVKfactor){
    zewkhists.push_back((TH2*) kffile_zewk->Get("TH2F_kFactor"));
    wewkhists.push_back((TH2*) kffile_wewk->Get("TH2F_kFactor"));
  }

  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;

  if(not nloSamples.useDYJetsNLO){
    cout<<"lepton+jets control region --> Z+jets"<<endl;
    makehist4(vlltree,vllhist,vllhist_2D,true,sample,category,false,1.0,lumi,vllhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }
  else{
    cout<<"lepton+jets control region --> Z+jets"<<endl;
    makehist4(vlltree_nlo1,vllhist,vllhist_2D,true,sample,category,false,1.00,lumi,vllhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree_nlo2,vllhist,vllhist_2D,true,sample,category,false,1.00,lumi,vllhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree_nlo3,vllhist,vllhist_2D,true,sample,category,false,1.00,lumi,vllhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree_nlo4,vllhist,vllhist_2D,true,sample,category,false,1.00,lumi,vllhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }

  if(not nloSamples.useWJetsNLO){
    cout<<"lepton+jets control region --> W+jets"<<endl;
    makehist4(vltree,vlhist,vlhist_2D,true,sample,category,false,1.0,lumi,vlhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }
  else{
    cout<<"lepton+jets control region --> W+jets"<<endl;
    makehist4(vltree_nlo1,vlhist,vlhist_2D,true,sample,category,false,1.00,lumi,vlhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree_nlo2,vlhist,vlhist_2D,true,sample,category,false,1.00,lumi,vlhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree_nlo3,vlhist,vlhist_2D,true,sample,category,false,1.00,lumi,vlhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree_nlo4,vlhist,vlhist_2D,true,sample,category,false,1.00,lumi,vlhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }

  cout<<"lepton+jets control region --> top"<<endl;
  makehist4(tttree,tthist,tthist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"lepton+jets control region --> Diboson"<<endl;
  makehist4(dbtree,dbhist,dbhist_2D,true,sample,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"lepton+jets control region --> gamma+jets"<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,sample,category,false,1.00,lumi,ahists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  cout<<"lepton+jets control region --> QCD"<<endl;
  makehist4(qctree,qchist,qchist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"lepton+jets control region --> EWK W"<<endl;
  makehist4(ewkwtree,ewkwhist,ewkwhist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight, wewkhists);
  cout<<"lepton+jets control region --> EWK Z"<<endl;
  makehist4(ewkztree,ewkzhist,ewkzhist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight, zewkhists);
  cout<<"lepton+jets control region --> ZllGamma"<<endl;
  makehist4(zgtree,vghist,vghist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"lepton+jets control region --> Wgamma"<<endl;
  makehist4(wgtree,vghist,vghist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  if(doShapeSystematics and (sample == Sample::zmm or sample == Sample::zee)){
    cout<<"lepton +jets region --> systematics for W+jets"<<endl;
    makehist4(vltree,vlhist_metJetUp,vlhist_metJetUp_2D,true,sample,category,false,1.00,lumi,vlhists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metJetDw,vlhist_metJetDw_2D,true,sample,category,false,1.00,lumi,vlhists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metResUp,vlhist_metResUp_2D,true,sample,category,false,1.00,lumi,vlhists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metResDw,vlhist_metResDw_2D,true,sample,category,false,1.00,lumi,vlhists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metUncUp,vlhist_metUncUp_2D,true,sample,category,false,1.00,lumi,vlhists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metUncDw,vlhist_metUncDw_2D,true,sample,category,false,1.00,lumi,vlhists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);
  }
  else if(doShapeSystematics and (sample == Sample::wmn or sample == Sample::wen)){
    cout<<"lepton +jets region --> systematics for Z+jets"<<endl;
    makehist4(vlltree,vllhist_metJetUp,vllhist_metJetUp_2D,true,sample,category,false,1.00,lumi,vllhists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metJetDw,vllhist_metJetDw_2D,true,sample,category,false,1.00,lumi,vllhists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metResUp,vllhist_metResUp_2D,true,sample,category,false,1.00,lumi,vllhists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metResDw,vllhist_metResDw_2D,true,sample,category,false,1.00,lumi,vllhists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metUncUp,vllhist_metUncUp_2D,true,sample,category,false,1.00,lumi,vllhists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metUncDw,vllhist_metUncDw_2D,true,sample,category,false,1.00,lumi,vllhists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);
  }

  if(doShapeSystematics){

    cout<<"lepton +jets region --> systematics for top"<<endl;
    makehist4(tttree,tthist_metJetUp,tthist_metJetUp_2D,true,sample,category,false,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metJetDw,tthist_metJetDw_2D,true,sample,category,false,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metResUp,tthist_metResUp_2D,true,sample,category,false,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metResDw,tthist_metResDw_2D,true,sample,category,false,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metUncUp,tthist_metUncUp_2D,true,sample,category,false,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(tttree,tthist_metUncDw,tthist_metUncDw_2D,true,sample,category,false,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"lepton +jets region --> systematics for di-boson"<<endl;
    makehist4(dbtree,dbhist_metJetUp,dbhist_metJetUp_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metJetDw,dbhist_metJetDw_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metResUp,dbhist_metResUp_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metResDw,dbhist_metResDw_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metUncUp,dbhist_metUncUp_2D,true,sample,category,isWJet,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metUncDw,dbhist_metUncDw_2D,true,sample,category,isWJet,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"lepton +jets region --> systematics for gamma+jets"<<endl;
    makehist4(gmtree,gmhist_metJetUp,gmhist_metJetUp_2D,true,sample,category,false,1.00,lumi,ahists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metJetDw,gmhist_metJetDw_2D,true,sample,category,false,1.00,lumi,ahists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metResUp,gmhist_metResUp_2D,true,sample,category,false,1.00,lumi,ahists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metResDw,gmhist_metResDw_2D,true,sample,category,false,1.00,lumi,ahists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metUncUp,gmhist_metUncUp_2D,true,sample,category,false,1.00,lumi,ahists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metUncDw,gmhist_metUncDw_2D,true,sample,category,false,1.00,lumi,ahists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);


  }

  cout<<"lepton+jets control region --> Data"<<endl;
  makehist4(dttree_merged,dthist,dthist_2D,false,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  //smooth
  if(doSmoothing){
    for(auto hist : tthist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qchist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vghist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    
    for(auto hist : tthist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qchist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vghist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    

    if(doShapeSystematics){
      
      for(auto hist : tthist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : dbhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : vlhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

      for(auto hist : vllhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : gmhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : tthist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : tthist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : dbhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : dbhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : vlhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vlhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : vllhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : vllhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      
      for(auto hist : gmhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
      for(auto hist : gmhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    }
  }
  
  //
  string dirName;
  if(sample == Sample::zmm)
    dirName = "ZM";
  else if(sample == Sample::wmn)
    dirName = "WM";
  else if(sample == Sample::zee)
    dirName = "ZE";
  else if(sample == Sample::wen)
    dirName = "WE";

  outfile->cd();
  if(not outfile->GetDirectory(dirName.c_str()))
    outfile->mkdir(dirName.c_str());
  outfile->cd(dirName.c_str());
  for(auto hist :  dthist) hist->Write();
  for(auto hist :  tthist) hist->Write();
  for(auto hist :  dbhist) hist->Write();
  for(auto hist :  gmhist) hist->Write();
  for(auto hist :  qchist) hist->Write();
  for(auto hist :  vlhist) hist->Write();
  for(auto hist :  vllhist) hist->Write();
  for(auto hist :  ewkwhist) hist->Write();
  for(auto hist :  ewkzhist) hist->Write();
  for(auto hist :  vghist) hist->Write();

  if(doShapeSystematics){
    outfile->cd();
    if(not outfile->GetDirectory((dirName+"/sysShape").c_str()))
      outfile->mkdir((dirName+"/sysShape").c_str());
    outfile->cd((dirName+"/sysShape").c_str());

    for(auto hist :  dbhist_metJetUp) hist->Write();
    for(auto hist :  dbhist_metJetDw) hist->Write();
    for(auto hist :  dbhist_metResUp) hist->Write();
    for(auto hist :  dbhist_metResDw) hist->Write();
    for(auto hist :  dbhist_metUncUp) hist->Write();
    for(auto hist :  dbhist_metUncDw) hist->Write();

    for(auto hist :  gmhist_metJetUp) hist->Write();
    for(auto hist :  gmhist_metJetDw) hist->Write();
    for(auto hist :  gmhist_metResUp) hist->Write();
    for(auto hist :  gmhist_metResDw) hist->Write();
    for(auto hist :  gmhist_metUncUp) hist->Write();
    for(auto hist :  gmhist_metUncDw) hist->Write();
 
    for(auto hist :  tthist_metJetUp) hist->Write();
    for(auto hist :  tthist_metJetDw) hist->Write();
    for(auto hist :  tthist_metResUp) hist->Write();
    for(auto hist :  tthist_metResDw) hist->Write();
    for(auto hist :  tthist_metUncUp) hist->Write();
    for(auto hist :  tthist_metUncDw) hist->Write();

    for(auto hist :  vlhist_metJetUp) hist->Write();
    for(auto hist :  vlhist_metJetDw) hist->Write();
    for(auto hist :  vlhist_metResUp) hist->Write();
    for(auto hist :  vlhist_metResDw) hist->Write();
    for(auto hist :  vlhist_metUncUp) hist->Write();
    for(auto hist :  vlhist_metUncDw) hist->Write();

    for(auto hist :  vllhist_metJetUp) hist->Write();
    for(auto hist :  vllhist_metJetDw) hist->Write();
    for(auto hist :  vllhist_metResUp) hist->Write();
    for(auto hist :  vllhist_metResDw) hist->Write();
    for(auto hist :  vllhist_metUncUp) hist->Write();
    for(auto hist :  vllhist_metUncDw) hist->Write();
  }

  outfile->cd();
  outfile->cd(dirName.c_str());

  for(auto hist_2D : vlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vllhist_2D){TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : tthist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dbhist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : gmhist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : qchist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dthist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkzhist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkwhist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vghist_2D) {TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

  if(doShapeSystematics){
    
    outfile->cd();
    outfile->cd((dirName+"/sysShape").c_str());
    for(auto hist : dbhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : dbhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

    for(auto hist : gmhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : gmhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
 
    for(auto hist : tthist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : tthist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

    for(auto hist : vlhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vlhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

    for(auto hist : vllhist_metJetUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metJetDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metResUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metResDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metUncUp_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }
    for(auto hist : vllhist_metUncDw_2D){ TH1* temp = unroll2DHistograms(hist); temp->Write(); }

  }
  
  outfile->cd();

  vlhists.clear();
  vllhists.clear();
  ahists.clear();
  zewkhists.clear();
  wewkhists.clear();

  for(auto hist : dthist) delete hist;
  for(auto hist : tthist) delete hist;
  for(auto hist : tthist_metJetUp) delete hist;
  for(auto hist : tthist_metJetDw) delete hist;
  for(auto hist : tthist_metResUp) delete hist;
  for(auto hist : tthist_metResDw) delete hist;
  for(auto hist : tthist_metUncUp) delete hist;
  for(auto hist : tthist_metUncDw) delete hist;
  for(auto hist : dbhist) delete hist;
  for(auto hist : dbhist_metJetUp) delete hist;
  for(auto hist : dbhist_metJetDw) delete hist;
  for(auto hist : dbhist_metResUp) delete hist;
  for(auto hist : dbhist_metResDw) delete hist;
  for(auto hist : dbhist_metUncUp) delete hist;
  for(auto hist : dbhist_metUncDw) delete hist;
  for(auto hist : qchist) delete hist;
  for(auto hist : gmhist) delete hist;
  for(auto hist : gmhist_metJetUp) delete hist;
  for(auto hist : gmhist_metJetDw) delete hist;
  for(auto hist : gmhist_metResUp) delete hist;
  for(auto hist : gmhist_metResDw) delete hist;
  for(auto hist : gmhist_metUncUp) delete hist;
  for(auto hist : gmhist_metUncDw) delete hist;  
  for(auto hist : vlhist) delete hist;
  for(auto hist : vlhist_metJetUp) delete hist;
  for(auto hist : vlhist_metJetDw) delete hist;
  for(auto hist : vlhist_metResUp) delete hist;
  for(auto hist : vlhist_metResDw) delete hist;
  for(auto hist : vlhist_metUncUp) delete hist;
  for(auto hist : vlhist_metUncDw) delete hist;
  for(auto hist : vllhist) delete hist;
  for(auto hist : vllhist_metJetUp) delete hist;
  for(auto hist : vllhist_metJetDw) delete hist;
  for(auto hist : vllhist_metResUp) delete hist;
  for(auto hist : vllhist_metResDw) delete hist;
  for(auto hist : vllhist_metUncUp) delete hist;
  for(auto hist : vllhist_metUncDw) delete hist;
  for(auto hist : ewkzhist) delete hist;
  for(auto hist : ewkwhist) delete hist;
  for(auto hist : vghist) delete hist;
  
  dthist.clear();
  tthist.clear();
  tthist_metJetUp.clear();
  tthist_metJetDw.clear();
  tthist_metResUp.clear();
  tthist_metResDw.clear();
  tthist_metUncUp.clear();
  tthist_metUncDw.clear();
  qchist.clear();
  dbhist.clear();
  dbhist_metJetUp.clear();
  dbhist_metJetDw.clear();
  dbhist_metResUp.clear();
  dbhist_metResDw.clear();
  dbhist_metUncUp.clear();
  dbhist_metUncDw.clear();
  gmhist.clear();
  gmhist_metJetUp.clear();
  gmhist_metJetDw.clear();
  gmhist_metResUp.clear();
  gmhist_metResDw.clear();
  gmhist_metUncUp.clear();
  gmhist_metUncDw.clear();
  vlhist.clear();
  vlhist_metJetUp.clear();
  vlhist_metJetDw.clear();
  vlhist_metResUp.clear();
  vlhist_metResDw.clear();
  vlhist_metUncUp.clear();
  vlhist_metUncDw.clear();
  vllhist.clear();
  vllhist_metJetUp.clear();
  vllhist_metJetDw.clear();
  vllhist_metResUp.clear();
  vllhist_metResDw.clear();
  vllhist_metUncUp.clear();
  vllhist_metUncDw.clear();
  ewkwhist.clear();
  ewkzhist.clear();
  vghist.clear();

  for(auto hist_2D : dthist_2D) delete hist_2D;
  for(auto hist_2D : tthist_2D) delete hist_2D;
  for(auto hist_2D : tthist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D : tthist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D : tthist_metResUp_2D) delete hist_2D;
  for(auto hist_2D : tthist_metResDw_2D) delete hist_2D;
  for(auto hist_2D : tthist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D : tthist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D : dbhist_2D) delete hist_2D;
  for(auto hist_2D : dbhist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D : dbhist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D : dbhist_metResUp_2D) delete hist_2D;
  for(auto hist_2D : dbhist_metResDw_2D) delete hist_2D;
  for(auto hist_2D : dbhist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D : dbhist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D : qchist_2D) delete hist_2D;
  for(auto hist_2D : gmhist_2D) delete hist_2D;
  for(auto hist_2D : gmhist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D : gmhist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D : gmhist_metResUp_2D) delete hist_2D;
  for(auto hist_2D : gmhist_metResDw_2D) delete hist_2D;
  for(auto hist_2D : gmhist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D : gmhist_metUncDw_2D) delete hist_2D;  
  for(auto hist_2D : vlhist_2D) delete hist_2D;
  for(auto hist_2D : vlhist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D : vlhist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D : vlhist_metResUp_2D) delete hist_2D;
  for(auto hist_2D : vlhist_metResDw_2D) delete hist_2D;
  for(auto hist_2D : vlhist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D : vlhist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D : vllhist_2D) delete hist_2D;
  for(auto hist_2D : vllhist_metJetUp_2D) delete hist_2D;
  for(auto hist_2D : vllhist_metJetDw_2D) delete hist_2D;
  for(auto hist_2D : vllhist_metResUp_2D) delete hist_2D;
  for(auto hist_2D : vllhist_metResDw_2D) delete hist_2D;
  for(auto hist_2D : vllhist_metUncUp_2D) delete hist_2D;
  for(auto hist_2D : vllhist_metUncDw_2D) delete hist_2D;
  for(auto hist_2D : ewkzhist_2D) delete hist_2D;
  for(auto hist_2D : ewkwhist_2D) delete hist_2D;
  for(auto hist_2D : vghist_2D) delete hist_2D;

  dthist_2D.clear();
  tthist_2D.clear();
  tthist_metJetUp_2D.clear();
  tthist_metJetDw_2D.clear();
  tthist_metResUp_2D.clear();
  tthist_metResDw_2D.clear();
  tthist_metUncUp_2D.clear();
  tthist_metUncDw_2D.clear();
  qchist_2D.clear();
  dbhist_2D.clear();
  dbhist_metJetUp_2D.clear();
  dbhist_metJetDw_2D.clear();
  dbhist_metResUp_2D.clear();
  dbhist_metResDw_2D.clear();
  dbhist_metUncUp_2D.clear();
  dbhist_metUncDw_2D.clear();
  gmhist_2D.clear();
  gmhist_metJetUp_2D.clear();
  gmhist_metJetDw_2D.clear();
  gmhist_metResUp_2D.clear();
  gmhist_metResDw_2D.clear();
  gmhist_metUncUp_2D.clear();
  gmhist_metUncDw_2D.clear();
  vlhist_2D.clear();
  vlhist_metJetUp_2D.clear();
  vlhist_metJetDw_2D.clear();
  vlhist_metResUp_2D.clear();
  vlhist_metResDw_2D.clear();
  vlhist_metUncUp_2D.clear();
  vlhist_metUncDw_2D.clear();
  vllhist_2D.clear();
  vllhist_metJetUp_2D.clear();
  vllhist_metJetDw_2D.clear();
  vllhist_metResUp_2D.clear();
  vllhist_metResDw_2D.clear();
  vllhist_metUncUp_2D.clear();
  vllhist_metUncDw_2D.clear();

  ewkwhist_2D.clear();
  ewkzhist_2D.clear();
  vghist_2D.clear();
  
  if(kfactzjet_vbf)
    kfactzjet_vbf->Close();
  if(kfactwjet_vbf)
    kfactwjet_vbf->Close();
  if(kfactgjet_vbf)
    kfactgjet_vbf->Close();
  if(kffile)
     kffile->Close();
  if(kffile_alt)
    kffile_alt->Close();
  if(kffile_unlops)
    kffile_unlops->Close();
  if(kffile_zewk)
    kffile_zewk->Close();
  if(kffile_wewk)
    kffile_wewk->Close();

  if(tttree) delete tttree;
  if(dbtree) delete dbtree;
  if(gmtree) delete gmtree;
  if(qctree) delete qctree;
  if(vltree) delete vltree;
  if(vlltree) delete vlltree;
  if(ewkwtree) delete ewkwtree;
  if(ewkztree) delete ewkztree;
  if(zgtree) delete zgtree;
  if(wgtree) delete wgtree;
  if(vltree_nlo1) delete vltree_nlo1;
  if(vltree_nlo2) delete vltree_nlo2;
  if(vltree_nlo3) delete vltree_nlo3;
  if(vltree_nlo4) delete vltree_nlo4;
  if(vlltree_nlo1) delete vlltree_nlo1;
  if(vlltree_nlo2) delete vlltree_nlo2;
  if(vlltree_nlo3) delete vlltree_nlo3;
  if(vlltree_nlo4) delete vlltree_nlo4;
  if(dttree) delete dttree;
  if(dttree_2) delete dttree_2;

  cout << "Templates for the lepton control region computed ..." << endl;
}

//build template for tt
void topdatamchist(TFile* outfile,
		   const Sample   & sample,
		   const Category & category,
		   vector<string> observables, vector<string> observables_2D,
		   const double & lumi,
		   const SamplesNLO & nloSamples,
		   const bool & makeResonantSelection = false,
		   const bool & doShapeSystematics    = false,
		   const bool & isHInv                = false,
		   const bool & doAlternativeTop      = false,
		   const bool & applyPFWeight         = false) {

  if (sample != Sample::topmu && sample != Sample::topel) return;

  TChain* tttree      = new TChain("tree/tree");
  TChain* tttree_alt  = new TChain("tree/tree");
  TChain* dbtree  = new TChain("tree/tree");
  TChain* gmtree  = new TChain("tree/tree");
  TChain* qctree  = new TChain("tree/tree");
  TChain* vltree  = new TChain("tree/tree");
  TChain* vlltree = new TChain("tree/tree");
  TChain* dttree  = new TChain("tree/tree");

  string suffix;

  if(sample == Sample::topmu)
    suffix = "topmu";
  else if(sample == Sample::topel)
    suffix = "topel";

  vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/"+suffix+"filter/*root").c_str());
  vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/"+suffix+"filter/*root").c_str());
  qctree->Add((baseInputTreePath+"/QCD/"+suffix+"filter/*root").c_str());
  dbtree->Add((baseInputTreePath+"/DiBoson/"+suffix+"filter/*root").c_str());
  gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/"+suffix+"filter/*root").c_str());
  tttree->Add((baseInputTreePath+"/TopAlternative/"+suffix+"filter/*root").c_str());
  tttree_alt->Add((baseInputTreePath+"/Top/"+suffix+"filter/*root").c_str());
  if(sample == Sample::topmu)
    //    dttree->Add((baseInputTreePath+"/MET/"+suffix+"filter/*root").c_str());
    dttree->Add((baseInputTreePath+"/MET_reMiniAOD/"+suffix+"filter/*root").c_str());
  else if(sample == Sample::topel)
    dttree->Add((baseInputTreePath+"/SingleElectron/"+suffix+"filter/*root").c_str());
  
  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_alt;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> dbhist_metJetUp;
  vector<TH1*> dbhist_metJetDw;
  vector<TH1*> dbhist_metResUp;
  vector<TH1*> dbhist_metResDw;
  vector<TH1*> dbhist_metUncUp;
  vector<TH1*> dbhist_metUncDw;

  vector<TH1*> gmhist;
  vector<TH1*> gmhist_metJetUp;
  vector<TH1*> gmhist_metJetDw;
  vector<TH1*> gmhist_metResUp;
  vector<TH1*> gmhist_metResDw;
  vector<TH1*> gmhist_metUncUp;
  vector<TH1*> gmhist_metUncDw;

  vector<TH1*> vlhist;
  vector<TH1*> vlhist_metJetUp;
  vector<TH1*> vlhist_metJetDw;
  vector<TH1*> vlhist_metResUp;
  vector<TH1*> vlhist_metResDw;
  vector<TH1*> vlhist_metUncUp;
  vector<TH1*> vlhist_metUncDw;
  vector<TH1*> vllhist;
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

  vector<double> bins;

  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable "+obs+"--> please define it"<<endl;

    TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* gmhist_temp = new TH1F((string("gmkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));

    if(makeResonantSelection){
      TH1F* tthist_matched_temp = new TH1F((string("tbkghist_matched")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_matched_alt_temp = new TH1F((string("tbkghist_matched_alt")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_unmatched_temp = new TH1F((string("tbkghist_unmatched")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_unmatched_alt_temp = new TH1F((string("tbkghist_unmatched_alt")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      tthist_matched.push_back(dynamic_cast<TH1*>(tthist_matched_temp));
      tthist_matched_alt.push_back(dynamic_cast<TH1*>(tthist_matched_alt_temp));
      tthist_unmatched.push_back(dynamic_cast<TH1*>(tthist_unmatched_temp));
      tthist_unmatched_alt.push_back(dynamic_cast<TH1*>(tthist_unmatched_alt_temp));

    }
    else{
      TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* tthist_alt_temp = new TH1F((string("tbkghist_alt")+suffix+"_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      tthist_alt.push_back(dynamic_cast<TH1*>(tthist_alt_temp));
    }

    if(doShapeSystematics){
      
      TH1F* dbhist_metJetUp_temp = new TH1F((string("dbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metJetDw_temp = new TH1F((string("dbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metResUp_temp = new TH1F((string("dbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metResDw_temp = new TH1F((string("dbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metUncUp_temp = new TH1F((string("dbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* dbhist_metUncDw_temp = new TH1F((string("dbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      
      dbhist_metJetUp.push_back(dynamic_cast<TH1*>(dbhist_metJetUp_temp));
      dbhist_metJetDw.push_back(dynamic_cast<TH1*>(dbhist_metJetDw_temp));
      dbhist_metResUp.push_back(dynamic_cast<TH1*>(dbhist_metResUp_temp));
      dbhist_metResDw.push_back(dynamic_cast<TH1*>(dbhist_metResDw_temp));
      dbhist_metUncUp.push_back(dynamic_cast<TH1*>(dbhist_metUncUp_temp));
      dbhist_metUncDw.push_back(dynamic_cast<TH1*>(dbhist_metUncDw_temp));


      TH1F* gmhist_metJetUp_temp = new TH1F((string("gbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metJetDw_temp = new TH1F((string("gbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metResUp_temp = new TH1F((string("gbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metResDw_temp = new TH1F((string("gbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metUncUp_temp = new TH1F((string("gbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* gmhist_metUncDw_temp = new TH1F((string("gbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      
      gmhist_metJetUp.push_back(dynamic_cast<TH1*>(gmhist_metJetUp_temp));
      gmhist_metJetDw.push_back(dynamic_cast<TH1*>(gmhist_metJetDw_temp));
      gmhist_metResUp.push_back(dynamic_cast<TH1*>(gmhist_metResUp_temp));
      gmhist_metResDw.push_back(dynamic_cast<TH1*>(gmhist_metResDw_temp));
      gmhist_metUncUp.push_back(dynamic_cast<TH1*>(gmhist_metUncUp_temp));
      gmhist_metUncDw.push_back(dynamic_cast<TH1*>(gmhist_metUncDw_temp));

      
      TH1F* vlhist_metJetUp_temp = new TH1F((string("vlbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metJetDw_temp = new TH1F((string("vlbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metResUp_temp = new TH1F((string("vlbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metResDw_temp = new TH1F((string("vlbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metUncUp_temp = new TH1F((string("vlbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vlhist_metUncDw_temp = new TH1F((string("vlbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

      vlhist_metJetUp.push_back(dynamic_cast<TH1*>(vlhist_metJetUp_temp));
      vlhist_metJetDw.push_back(dynamic_cast<TH1*>(vlhist_metJetDw_temp));
      vlhist_metResUp.push_back(dynamic_cast<TH1*>(vlhist_metResUp_temp));
      vlhist_metResDw.push_back(dynamic_cast<TH1*>(vlhist_metResDw_temp));
      vlhist_metUncUp.push_back(dynamic_cast<TH1*>(vlhist_metUncUp_temp));
      vlhist_metUncDw.push_back(dynamic_cast<TH1*>(vlhist_metUncDw_temp));

      TH1F* vllhist_metJetUp_temp = new TH1F((string("vllbkghist")+suffix+"_metJetUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metJetDw_temp = new TH1F((string("vllbkghist")+suffix+"_metJetDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metResUp_temp = new TH1F((string("vllbkghist")+suffix+"_metResUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metResDw_temp = new TH1F((string("vllbkghist")+suffix+"_metResDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metUncUp_temp = new TH1F((string("vllbkghist")+suffix+"_metUncUp_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      TH1F* vllhist_metUncDw_temp = new TH1F((string("vllbkghist")+suffix+"_metUncDw_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
 
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
  vector<TH2*> dbhist_metJetUp_2D;
  vector<TH2*> dbhist_metJetDw_2D;
  vector<TH2*> dbhist_metResUp_2D;
  vector<TH2*> dbhist_metResDw_2D;
  vector<TH2*> dbhist_metUncUp_2D;
  vector<TH2*> dbhist_metUncDw_2D;
  vector<TH2*> gmhist_2D;
  vector<TH2*> gmhist_metJetUp_2D;
  vector<TH2*> gmhist_metJetDw_2D;
  vector<TH2*> gmhist_metResUp_2D;
  vector<TH2*> gmhist_metResDw_2D;
  vector<TH2*> gmhist_metUncUp_2D;
  vector<TH2*> gmhist_metUncDw_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vlhist_metJetUp_2D;
  vector<TH2*> vlhist_metJetDw_2D;
  vector<TH2*> vlhist_metResUp_2D;
  vector<TH2*> vlhist_metResDw_2D;
  vector<TH2*> vlhist_metUncUp_2D;
  vector<TH2*> vlhist_metUncDw_2D;
  vector<TH2*> vllhist_2D;
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
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;

    TH2F* dthist_temp     = new TH2F((string("datahist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dbhist_temp     = new TH2F((string("dbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* gmhist_temp     = new TH2F((string("gbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* qchist_temp     = new TH2F((string("qbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vlhist_temp     = new TH2F((string("vlbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vllhist_temp    = new TH2F((string("vllbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* tthist_temp     = new TH2F((string("tbkghist")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* tthist_alt_temp = new TH2F((string("tbkghist_alt")+suffix+"_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
  
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    tthist_alt_2D.push_back(dynamic_cast<TH2*>(tthist_alt_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));
    dbhist_2D.push_back(dynamic_cast<TH2*>(dbhist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    qchist_2D.push_back(dynamic_cast<TH2*>(qchist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vllhist_2D.push_back(dynamic_cast<TH2*>(vllhist_temp));

    if(doShapeSystematics){
      
      TH2F* dbhist_metJetUp_temp = new TH2F((string("dbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metJetDw_temp = new TH2F((string("dbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metResUp_temp = new TH2F((string("dbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metResDw_temp = new TH2F((string("dbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metUncUp_temp = new TH2F((string("dbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* dbhist_metUncDw_temp = new TH2F((string("dbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      
      dbhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetUp_temp));
      dbhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metJetDw_temp));
      dbhist_metResUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metResUp_temp));
      dbhist_metResDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metResDw_temp));
      dbhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncUp_temp));
      dbhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(dbhist_metUncDw_temp));

      TH2F* gmhist_metJetUp_temp = new TH2F((string("gbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metJetDw_temp = new TH2F((string("gbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metResUp_temp = new TH2F((string("gbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metResDw_temp = new TH2F((string("gbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metUncUp_temp = new TH2F((string("gbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* gmhist_metUncDw_temp = new TH2F((string("gbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      
      gmhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetUp_temp));
      gmhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metJetDw_temp));
      gmhist_metResUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metResUp_temp));
      gmhist_metResDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metResDw_temp));
      gmhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncUp_temp));
      gmhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(gmhist_metUncDw_temp));
      
      TH2F* vlhist_metJetUp_temp = new TH2F((string("vlbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metJetDw_temp = new TH2F((string("vlbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metResUp_temp = new TH2F((string("vlbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metResDw_temp = new TH2F((string("vlbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metUncUp_temp = new TH2F((string("vlbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vlhist_metUncDw_temp = new TH2F((string("vlbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

      vlhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetUp_temp));
      vlhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metJetDw_temp));
      vlhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metResUp_temp));
      vlhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metResDw_temp));
      vlhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncUp_temp));
      vlhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vlhist_metUncDw_temp));

      TH2F* vllhist_metJetUp_temp = new TH2F((string("vllbkghist")+suffix+"_metJetUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metJetDw_temp = new TH2F((string("vllbkghist")+suffix+"_metJetDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metResUp_temp = new TH2F((string("vllbkghist")+suffix+"_metResUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metResDw_temp = new TH2F((string("vllbkghist")+suffix+"_metResDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metUncUp_temp = new TH2F((string("vllbkghist")+suffix+"_metUncUp_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      TH2F* vllhist_metUncDw_temp = new TH2F((string("vllbkghist")+suffix+"_metUncDw_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
 
      vllhist_metJetUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetUp_temp));
      vllhist_metJetDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metJetDw_temp));
      vllhist_metResUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metResUp_temp));
      vllhist_metResDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metResDw_temp));
      vllhist_metUncUp_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncUp_temp));
      vllhist_metUncDw_2D.push_back(dynamic_cast<TH2*>(vllhist_metUncDw_temp));
      
    }
  }

  // k-factors
  TFile kffile(kfactorFile.c_str());  
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1*  zewkhist  = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  TH1* wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1* wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1* wewkhist  = (TH1*) kffile.Get("EWKcorr/W");

  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1* alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  vector<TH1*> ehists;
  vector<TH1*> vlhists;
  vector<TH1*> ahists;
  vector<TH1*> vllhists;
  // apply NLO QCD and EWK corrections for Zll and Wlnu                                                                                                                      
  if(nloSamples.usePhotonJetsNLO)
    ahists.push_back(aewkhist);
  else{
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);
  }

  if(nloSamples.useWJetsNLO)
    vlhists.push_back(wewkhist);
  else{
    vlhists.push_back(wnlohist);
    vlhists.push_back(wewkhist);
  }
  
  if(nloSamples.useDYJetsNLO)
    vllhists.push_back(zewkhist);
  else{
    vllhists.push_back(znlohist);
    vllhists.push_back(zewkhist);
  }

  // 
  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;

  bool topPtReweight = false;
  if(category == Category::monoV or category == Category::boosted)
    topPtReweight = true;


  if(not makeResonantSelection){
    cout<<"top+jets control region --> Top"<<endl;
    makehist4(tttree,tthist,tthist_2D,true,sample,category,false,1.00,lumi,ehists,"",topPtReweight,reweightNVTX,0,isHInv,applyPFWeight);
    cout<<"top+jets control region --> Top alternative"<<endl;
    if(doAlternativeTop){
      makehist4(tttree_alt,tthist_alt,tthist_alt_2D,true,sample,category,false,1.00,lumi,ehists,"",topPtReweight,reweightNVTX,0,isHInv,applyPFWeight);
    }
  }
  else{
    cout<<"top+jets control region --> Top"<<endl;
    makehist4(tttree,tthist_matched,tthist_matched_2D,true,sample,category,false,1.00,lumi,ehists,"",topPtReweight,reweightNVTX,1,isHInv,applyPFWeight);
    makehist4(tttree,tthist_unmatched,tthist_unmatched_2D,true,sample,category,false,1.00,lumi,ehists,"",topPtReweight,reweightNVTX,2,isHInv,applyPFWeight);
    cout<<"top+jets control region --> Top alternative"<<endl;
    if(doAlternativeTop){
      makehist4(tttree_alt,tthist_matched_alt,tthist_matched_alt_2D,true,sample,category,false,1.00,lumi,ehists,"",topPtReweight,reweightNVTX,1,isHInv,applyPFWeight);
      makehist4(tttree_alt,tthist_unmatched_alt,tthist_unmatched_alt_2D,true,sample,category,false,1.00,lumi,ehists,"",topPtReweight,reweightNVTX,2,isHInv,applyPFWeight);
    }
  }

  cout<<"top+jets control region --> W+jets"<<endl;
  makehist4(vltree,vlhist,vlhist_2D,true,sample,category,false,1.00,lumi,vlhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"top+jets control region --> Z+jets"<<endl;
  makehist4(vlltree,vllhist,vllhist_2D,true,sample,category,false,1.00,lumi,vllhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"top+jets control region: Diboson"<<endl;
  makehist4(dbtree,dbhist,dbhist_2D,true,sample,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"top+jets control region: gamma+jets"<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,sample,category,false,1.00,lumi,ahists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"top+jets control region: QCD"<<endl;
  makehist4(qctree,qchist,qchist_2D,true,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  if(doShapeSystematics){
    cout<<"top+jets control region -->  shape systematics for W+Jets"<<endl;
    makehist4(vltree,vlhist_metJetUp,vlhist_metJetUp_2D,true,sample,category,false,1.00,lumi,vlhists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metJetDw,vlhist_metJetDw_2D,true,sample,category,false,1.00,lumi,vlhists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metResUp,vlhist_metResUp_2D,true,sample,category,false,1.00,lumi,vlhists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metResDw,vlhist_metResDw_2D,true,sample,category,false,1.00,lumi,vlhists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metUncUp,vlhist_metUncUp_2D,true,sample,category,false,1.00,lumi,vlhists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vltree,vlhist_metUncDw,vlhist_metUncDw_2D,true,sample,category,false,1.00,lumi,vlhists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"top+jets control region -->  shape systematics for Z+Jets"<<endl;
    makehist4(vlltree,vllhist_metJetUp,vllhist_metJetUp_2D,true,sample,category,false,1.00,lumi,vllhists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metJetDw,vllhist_metJetDw_2D,true,sample,category,false,1.00,lumi,vllhists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metResUp,vllhist_metResUp_2D,true,sample,category,false,1.00,lumi,vllhists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metResDw,vllhist_metResDw_2D,true,sample,category,false,1.00,lumi,vllhists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metUncUp,vllhist_metUncUp_2D,true,sample,category,false,1.00,lumi,vllhists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(vlltree,vllhist_metUncDw,vllhist_metUncDw_2D,true,sample,category,false,1.00,lumi,vllhists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"top+jets control region -->  shape systematics for Dibosons"<<endl;
    makehist4(dbtree,dbhist_metJetUp,dbhist_metJetUp_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metJetDw,dbhist_metJetDw_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metResUp,dbhist_metResUp_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metResDw,dbhist_metResDw_2D,true,sample,category,isWJet,1.00,lumi,ehists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metUncUp,dbhist_metUncUp_2D,true,sample,category,isWJet,1.00,lumi,ehists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(dbtree,dbhist_metUncDw,dbhist_metUncDw_2D,true,sample,category,isWJet,1.00,lumi,ehists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);

    cout<<"top+jets control region -->  shape systematics for gamma+jets"<<endl;
    makehist4(gmtree,gmhist_metJetUp,gmhist_metJetUp_2D,true,sample,category,false,1.00,lumi,ahists,"jesUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metJetDw,gmhist_metJetDw_2D,true,sample,category,false,1.00,lumi,ahists,"jesDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metResUp,gmhist_metResUp_2D,true,sample,category,false,1.00,lumi,ahists,"jerUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metResDw,gmhist_metResDw_2D,true,sample,category,false,1.00,lumi,ahists,"jerDw",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metUncUp,gmhist_metUncUp_2D,true,sample,category,false,1.00,lumi,ahists,"uncUp",false,reweightNVTX,0,isHInv,applyPFWeight);
    makehist4(gmtree,gmhist_metUncDw,gmhist_metUncDw_2D,true,sample,category,false,1.00,lumi,ahists,"uncDw",false,reweightNVTX,0,isHInv,applyPFWeight);
  }


  cout<<"top+jets control region: Data"<<endl;
  makehist4(dttree,dthist,dthist_2D,false,sample,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  
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

    for(auto hist : dbhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vlhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vllhist_metJetUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncUp){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncDw){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : dbhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : gmhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vlhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vlhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}

    for(auto hist : vllhist_metJetUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metJetDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metResDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncUp_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_metUncDw_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  }

  string dirName;
  if(sample == Sample::topmu)
    dirName = "TM";
  else if(sample == Sample::topel)
    dirName = "TE";

  outfile->cd();
  if(not outfile->GetDirectory((dirName).c_str()))
    outfile->mkdir(dirName.c_str());
  outfile->cd(dirName.c_str());

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
    outfile->cd();
    if(not outfile->GetDirectory((dirName+"/sysShape").c_str()))
      outfile->mkdir((dirName+"/sysShape").c_str());
    outfile->cd((dirName+"/sysShape").c_str());

    for(auto hist :  dbhist_metJetUp) hist->Write();
    for(auto hist :  dbhist_metJetDw) hist->Write();
    for(auto hist :  dbhist_metResUp) hist->Write();
    for(auto hist :  dbhist_metResDw) hist->Write();
    for(auto hist :  dbhist_metUncUp) hist->Write();
    for(auto hist :  dbhist_metUncDw) hist->Write();

    for(auto hist :  gmhist_metJetUp) hist->Write();
    for(auto hist :  gmhist_metJetDw) hist->Write();
    for(auto hist :  gmhist_metResUp) hist->Write();
    for(auto hist :  gmhist_metResDw) hist->Write();
    for(auto hist :  gmhist_metUncUp) hist->Write();
    for(auto hist :  gmhist_metUncDw) hist->Write();

    for(auto hist :  vlhist_metJetUp) hist->Write();
    for(auto hist :  vlhist_metJetDw) hist->Write();
    for(auto hist :  vlhist_metResUp) hist->Write();
    for(auto hist :  vlhist_metResDw) hist->Write();
    for(auto hist :  vlhist_metUncUp) hist->Write();
    for(auto hist :  vlhist_metUncDw) hist->Write();

    for(auto hist :  vllhist_metJetUp) hist->Write();
    for(auto hist :  vllhist_metJetDw) hist->Write();
    for(auto hist :  vllhist_metResUp) hist->Write();
    for(auto hist :  vllhist_metResDw) hist->Write();
    for(auto hist :  vllhist_metUncUp) hist->Write();
    for(auto hist :  vllhist_metUncDw) hist->Write();
  }
  
  outfile->cd();
  outfile->cd(dirName.c_str());
  for(auto hist_2D : vlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : vllhist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : tthist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : dbhist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : gmhist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : qchist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}
  for(auto hist_2D : dthist_2D) { TH1* temp = unroll2DHistograms(hist_2D); temp->Write();}

  if(doShapeSystematics){
    outfile->cd();
    outfile->cd((dirName+"/sysShape").c_str());

    for(auto hist :  dbhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  dbhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist :  gmhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  gmhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist :  vlhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vlhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}

    for(auto hist :  vllhist_metJetUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metJetDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metResUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metResDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metUncUp_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
    for(auto hist :  vllhist_metUncDw_2D) { TH1* temp = unroll2DHistograms(hist); temp->Write();}
  }

  outfile->cd();

  for(auto hist : dthist) delete hist;
  for(auto hist : tthist) delete hist;
  for(auto hist : tthist_alt) delete hist;
  for(auto hist : qchist) delete hist;
  for(auto hist : dbhist) delete hist;
  for(auto hist : dbhist_metJetUp) delete hist;
  for(auto hist : dbhist_metJetDw) delete hist;
  for(auto hist : dbhist_metResUp) delete hist;
  for(auto hist : dbhist_metResDw) delete hist;
  for(auto hist : dbhist_metUncUp) delete hist;
  for(auto hist : dbhist_metUncDw) delete hist;
  for(auto hist : gmhist) delete hist;
  for(auto hist : gmhist_metJetUp) delete hist;
  for(auto hist : gmhist_metJetDw) delete hist;
  for(auto hist : gmhist_metResUp) delete hist;
  for(auto hist : gmhist_metResDw) delete hist;
  for(auto hist : gmhist_metUncUp) delete hist;
  for(auto hist : gmhist_metUncDw) delete hist;
  for(auto hist : vlhist) delete hist;
  for(auto hist : vlhist_metJetUp) delete hist;
  for(auto hist : vlhist_metJetDw) delete hist;
  for(auto hist : vlhist_metResUp) delete hist;
  for(auto hist : vlhist_metResDw) delete hist;
  for(auto hist : vlhist_metUncUp) delete hist;
  for(auto hist : vlhist_metUncDw) delete hist;
  for(auto hist : vllhist) delete hist;
  for(auto hist : vllhist_metJetUp) delete hist;
  for(auto hist : vllhist_metJetDw) delete hist;
  for(auto hist : vllhist_metResUp) delete hist;
  for(auto hist : vllhist_metResDw) delete hist;
  for(auto hist : vllhist_metUncUp) delete hist;
  for(auto hist : vllhist_metUncDw) delete hist;

  dthist.clear();
  tthist.clear();
  tthist_alt.clear();
  qchist.clear();
  dbhist.clear();
  dbhist_metJetUp.clear();
  dbhist_metJetDw.clear();
  dbhist_metResUp.clear();
  dbhist_metResDw.clear();
  dbhist_metUncUp.clear();
  dbhist_metUncDw.clear();
  gmhist.clear();
  gmhist_metJetUp.clear();
  gmhist_metJetDw.clear();
  gmhist_metResUp.clear();
  gmhist_metResDw.clear();
  gmhist_metUncUp.clear();
  gmhist_metUncDw.clear();
  vlhist.clear();
  vlhist_metJetUp.clear();
  vlhist_metJetDw.clear();
  vlhist_metResUp.clear();
  vlhist_metResDw.clear();
  vlhist_metUncUp.clear();
  vlhist_metUncDw.clear();
  vllhist.clear();
  vllhist_metJetUp.clear();
  vllhist_metJetDw.clear();
  vllhist_metResUp.clear();
  vllhist_metResDw.clear();
  vllhist_metUncUp.clear();
  vllhist_metUncDw.clear();

  for(auto hist : tthist_matched) delete hist;
  for(auto hist : tthist_matched_alt) delete hist;
  for(auto hist : tthist_unmatched) delete hist;
  for(auto hist : tthist_unmatched_alt) delete hist;

  tthist_matched.clear();
  tthist_matched_alt.clear();
  tthist_unmatched.clear();
  tthist_unmatched_alt.clear();

  for(auto hist : dthist_2D) delete hist;
  for(auto hist : tthist_2D) delete hist;
  for(auto hist : tthist_alt_2D) delete hist;
  for(auto hist : qchist_2D) delete hist;
  for(auto hist : dbhist_2D) delete hist;
  for(auto hist : dbhist_metJetUp_2D) delete hist;
  for(auto hist : dbhist_metJetDw_2D) delete hist;
  for(auto hist : dbhist_metResUp_2D) delete hist;
  for(auto hist : dbhist_metResDw_2D) delete hist;
  for(auto hist : dbhist_metUncUp_2D) delete hist;
  for(auto hist : dbhist_metUncDw_2D) delete hist;
  for(auto hist : gmhist_2D) delete hist;
  for(auto hist : gmhist_metJetUp_2D) delete hist;
  for(auto hist : gmhist_metJetDw_2D) delete hist;
  for(auto hist : gmhist_metResUp_2D) delete hist;
  for(auto hist : gmhist_metResDw_2D) delete hist;
  for(auto hist : gmhist_metUncUp_2D) delete hist;
  for(auto hist : gmhist_metUncDw_2D) delete hist;
  for(auto hist : vlhist_2D) delete hist;
  for(auto hist : vlhist_metJetUp_2D) delete hist;
  for(auto hist : vlhist_metJetDw_2D) delete hist;
  for(auto hist : vlhist_metResUp_2D) delete hist;
  for(auto hist : vlhist_metResDw_2D) delete hist;
  for(auto hist : vlhist_metUncUp_2D) delete hist;
  for(auto hist : vlhist_metUncDw_2D) delete hist;
  for(auto hist : vllhist_2D) delete hist;
  for(auto hist : vllhist_metJetUp_2D) delete hist;
  for(auto hist : vllhist_metJetDw_2D) delete hist;
  for(auto hist : vllhist_metResUp_2D) delete hist;
  for(auto hist : vllhist_metResDw_2D) delete hist;
  for(auto hist : vllhist_metUncUp_2D) delete hist;
  for(auto hist : vllhist_metUncDw_2D) delete hist;


  dthist_2D.clear();
  tthist_2D.clear();
  tthist_alt_2D.clear();
  qchist_2D.clear();
  dbhist_2D.clear();
  dbhist_metJetUp_2D.clear();
  dbhist_metJetDw_2D.clear();
  dbhist_metResUp_2D.clear();
  dbhist_metResDw_2D.clear();
  dbhist_metUncUp_2D.clear();
  dbhist_metUncDw_2D.clear();

  gmhist_2D.clear();
  gmhist_metJetUp_2D.clear();
  gmhist_metJetDw_2D.clear();
  gmhist_metResUp_2D.clear();
  gmhist_metResDw_2D.clear();
  gmhist_metUncUp_2D.clear();
  gmhist_metUncDw_2D.clear();

  vlhist_2D.clear();
  vlhist_metJetUp_2D.clear();
  vlhist_metJetDw_2D.clear();
  vlhist_metResUp_2D.clear();
  vlhist_metResDw_2D.clear();
  vlhist_metUncUp_2D.clear();
  vlhist_metUncDw_2D.clear();
  vllhist_2D.clear();
  vllhist_metJetUp_2D.clear();
  vllhist_metJetDw_2D.clear();
  vllhist_metResUp_2D.clear();
  vllhist_metResDw_2D.clear();
  vllhist_metUncUp_2D.clear();
  vllhist_metUncDw_2D.clear();

  if(tttree) delete tttree;
  if(tttree_alt) delete tttree_alt;
  if(dbtree) delete dbtree;
  if(gmtree) delete gmtree;
  if(qctree) delete qctree;
  if(vltree) delete vltree;
  if(vlltree) delete vlltree;
  if(dttree) delete dttree;
  
  cout << "Templates for the top control region computed ..." << endl;
}

/////////
void qcddatamchist(TFile* outfile,
                   const Category & category,
                   vector<string> observables,
		   vector<string> observables_2D,
                   const double & lumi,
		   const SamplesNLO & nloSamples,
		   const bool & doAlternativeTop    = false,
		   const bool & isHInv              = false,
		   const bool & applyPFWeight       = false) {


  // Files for Znunu,Wlnu, Zll, top, qcd , diboson, signal, data                                                                                                            
  TChain* vnntree = new TChain("tree/tree");
  TChain* vltree  = new TChain("tree/tree");
  TChain* vlltree = new TChain("tree/tree");
  TChain* tttree  = new TChain("tree/tree");
  TChain* tttree_alt = NULL;
  if(doAlternativeTop)
    tttree_alt    = new TChain("tree/tree");
  TChain* dbtree  = new TChain("tree/tree");
  TChain* gmtree  = new TChain("tree/tree");
  TChain* qcdtree = new TChain("tree/tree");
  TChain* ewkwtree = new TChain("tree/tree");
  TChain* ewkztree = new TChain("tree/tree");
  TChain* dttree   = new TChain("tree/tree");

  vnntree->Add((baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/*root").c_str());
  vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/*root").c_str());
  vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/sigfilter/*root").c_str());
  tttree->Add((baseInputTreePath+"/TopAlternative/sigfilter/*root").c_str());
  if(doAlternativeTop and tttree_alt != NULL){
    tttree_alt->Add((baseInputTreePath+"Top/sigfilter/*root").c_str());
  }

  qcdtree->Add((baseInputTreePath+"/QCD/sigfilter/*root").c_str());
  dbtree->Add((baseInputTreePath+"/DiBoson/sigfilter/*root").c_str());
  gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/sigfilter/*root").c_str());
  ewkztree->Add((baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/*root").c_str());
  ewkwtree->Add((baseInputTreePath+"/WJetsEWK/sigfilter/*root").c_str());
  dttree->Add((baseInputTreePath+"/MET_reMiniAOD/sigfilter/*root").c_str());

  // make met histograms                                                                                                                                                        
  vector<TH1*> vnnhist;
  vector<TH1*> vlhist;
  vector<TH1*> vllhist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_alt;
  vector<TH1*> dbhist;
  vector<TH1*> gmhist;
  vector<TH1*> qcdhist;
  vector<TH1*> ewkwhist;
  vector<TH1*> ewkzhist;
  vector<TH1*> dthist;

  /////////
  vector<double> bins;
  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;

    TH1F* vnnhist_temp = new TH1F(("vnnbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vlhist_temp = new TH1F(("vlbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vllhist_temp = new TH1F(("vllbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* tthist_temp = new TH1F(("tbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dbhist_temp = new TH1F(("dbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* gmhist_temp = new TH1F(("gbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* qcdhist_temp = new TH1F(("qbkghistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkwhist_temp = new TH1F(("ewkbkgwhistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkzhist_temp = new TH1F(("ewkbkgzhistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dthist_temp = new TH1F(("datahistqcd_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

    vnnhist.push_back(dynamic_cast<TH1*>(vnnhist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    ewkwhist.push_back(dynamic_cast<TH1*>(ewkwhist_temp));
    ewkzhist.push_back(dynamic_cast<TH1*>(ewkzhist_temp));
    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    
    if(doAlternativeTop){
      TH1F* tthist_alt_temp = new TH1F(("tbkghistqcd_alt_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      tthist_alt.push_back(dynamic_cast<TH1*>(tthist_alt_temp));
    }    
  }

  vector<TH2*> vnnhist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vllhist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_alt_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> gmhist_2D;
  vector<TH2*> qcdhist_2D;
  vector<TH2*> ewkwhist_2D;
  vector<TH2*> ewkzhist_2D;
  vector<TH2*> dthist_2D;


  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable "+obs+"--> please define it"<<endl;
    
    TH2F* vnnhist_temp = new TH2F(("vnnbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vlhist_temp = new TH2F(("vlbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vllhist_temp = new TH2F(("vllbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* tthist_temp = new TH2F(("tbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dbhist_temp = new TH2F(("dbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* gmhist_temp = new TH2F(("gbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* qcdhist_temp = new TH2F(("qbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkwhist_temp = new TH2F(("ewkwbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkzhist_temp = new TH2F(("ewkzbkghistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dthist_temp = new TH2F(("datahistqcd_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

    vnnhist_2D.push_back(dynamic_cast<TH2*>(vnnhist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vllhist_2D.push_back(dynamic_cast<TH2*>(vllhist_temp));
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    qcdhist_2D.push_back(dynamic_cast<TH2*>(qcdhist_temp));
    dbhist_2D.push_back(dynamic_cast<TH2*>(dbhist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    ewkwhist_2D.push_back(dynamic_cast<TH2*>(ewkwhist_temp));
    ewkzhist_2D.push_back(dynamic_cast<TH2*>(ewkzhist_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));

    if(doAlternativeTop){
      TH2F* tthist_alt_temp = new TH2F(("tbkghist_alt_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      tthist_alt_2D.push_back(dynamic_cast<TH2*>(tthist_alt_temp));
    }    
  }

  
  // get k-factors NLO                                                                                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1*  zewkhist = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1*  wewkhist = (TH1*) kffile.Get("EWKcorr/W");
 
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> dyhists;
  vector<TH1*> whists;
  vector<TH1*> ahists;

  // apply EWK and QCD corrections                                                                                                                                              
  float scale = 1;
  if(nloSamples.useZJetsNLO){
    zhists.push_back(zewkhist);
    scale = 3;
  }
  else{
    zhists.push_back(znlohist); 
    zhists.push_back(zewkhist);
  }
  
  if(nloSamples.useDYJetsNLO)
    dyhists.push_back(zewkhist);
  else{
    dyhists.push_back(znlohist); 
    dyhists.push_back(zewkhist);
  }
  
  if(nloSamples.useWJetsNLO)
    whists.push_back(wewkhist);    
  else{
    whists.push_back(wnlohist); 
    whists.push_back(wewkhist);    
  }

  if(nloSamples.usePhotonJetsNLO)
    ahists.push_back(aewkhist);
  else{
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);
  }


  //////////////////
  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;
 
  cout<<"QCD region: Z->nunu sample "<<endl;
  if(not nloSamples.useZJetsNLO)
    makehist4(vnntree,vnnhist,vnnhist_2D,true,Sample::qcd,category,false,1.00,lumi,zhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  else
    makehist4(vnntree,vnnhist,vnnhist_2D,true,Sample::qcd,category,false,scale,lumi,zhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  cout<<"QCD region: W+jets sample "<<endl;
  makehist4(vltree,vlhist,vlhist_2D,true,Sample::qcd,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"QCD region: Z+jets sample "<<endl;
  makehist4(vlltree,vllhist,vllhist_2D,true,Sample::qcd,category,false,1.00,lumi,dyhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"QCD region: gamma+jets sample "<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,Sample::qcd,category,false,1.00,lumi,ahists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"QCD region: ewkw+jets sample "<<endl;
  makehist4(ewkwtree,ewkwhist,ewkwhist_2D,true,Sample::qcd,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"QCD region: ewkz+jets sample "<<endl;
  makehist4(ewkztree,ewkzhist,ewkzhist_2D,true,Sample::qcd,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight); // temp fix for a wrong xsec
  cout<<"QCD region: TTbar sample "<<endl;
  makehist4(tttree,tthist,tthist_2D,true,Sample::qcd,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

    //alternative ttbar             
  if(doAlternativeTop){
    cout<<"QCD region: TTbar alternative sample "<<endl;
    makehist4(tttree_alt,tthist_alt,tthist_alt_2D,true,Sample::qcd,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }

  cout<<"QCD region: Diboson sample "<<endl;
  makehist4(dbtree,dbhist,dbhist_2D,true,Sample::qcd,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"QCD region: QCD sample "<<endl;
  makehist4(qcdtree,qcdhist,qcdhist_2D,true,Sample::qcd,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  // take average of ttbar                                                                                                                                                     
  if(doAlternativeTop){
    for(size_t iHisto = 0; iHisto < tthist.size(); iHisto++)
      makeAverage(tthist.at(iHisto),tthist_alt.at(iHisto));      
    for(size_t iHisto = 0; iHisto < tthist_2D.size(); iHisto++)
      makeAverage(tthist_2D.at(iHisto),tthist_alt_2D.at(iHisto));
  }

  //smooth
  if(doSmoothing){
    for(auto hist : tthist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qcdhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    
    for(auto hist : tthist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qcdhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  }

  // data                                                
  cout<<"QCD region analysis --> loop on data "<<endl;
  makehist4(dttree,dthist,dthist_2D,false,Sample::qcd,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  
  outfile->cd();
  if(not outfile->GetDirectory("QCD"))
    outfile->mkdir("QCD");
  outfile->cd("QCD");
  // store histograms                                                                                                                                                        
  for(auto hist : vnnhist) hist->Write();
  for(auto hist : vlhist) hist->Write();
  for(auto hist : vllhist) hist->Write();
  for(auto hist : tthist) hist->Write();
  for(auto hist : dbhist) hist->Write();
  for(auto hist : gmhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();
  for(auto hist : ewkwhist) hist->Write();
  for(auto hist : ewkzhist) hist->Write();
  for(auto hist : dthist) hist->Write();

  // store hist_2Dograms                                                                                                                                                     
  outfile->cd();
  outfile->cd("QCD");
  for(auto hist_2D : vnnhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vllhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : gmhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : tthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dbhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : qcdhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkwhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkzhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

  outfile->cd();

  for(auto hist : vnnhist) delete hist;
  for(auto hist : vlhist) delete hist;
  for(auto hist : vllhist) delete hist;
  for(auto hist : tthist) delete hist;
  for(auto hist : dbhist) delete hist;
  for(auto hist : gmhist) delete hist;
  for(auto hist : qcdhist) delete hist;
  for(auto hist : ewkwhist) delete hist;
  for(auto hist : ewkzhist) delete hist;
  for(auto hist : dthist) delete hist;

  vnnhist.clear();
  vlhist.clear();
  vllhist.clear();
  tthist.clear();
  dbhist.clear();
  gmhist.clear();
  qcdhist.clear();
  ewkwhist.clear();
  ewkzhist.clear();
  dthist.clear();

  for(auto hist_2D : vnnhist_2D) delete hist_2D;
  for(auto hist_2D : vlhist_2D) delete hist_2D;
  for(auto hist_2D : vllhist_2D) delete hist_2D;
  for(auto hist_2D : tthist_2D) delete hist_2D;
  for(auto hist_2D : dbhist_2D) delete hist_2D;
  for(auto hist_2D : gmhist_2D) delete hist_2D;
  for(auto hist_2D : qcdhist_2D) delete hist_2D;
  for(auto hist_2D : ewkwhist_2D) delete hist_2D;
  for(auto hist_2D : ewkzhist_2D) delete hist_2D;
  for(auto hist_2D : dthist_2D) delete hist_2D;

  vnnhist_2D.clear();
  vlhist_2D.clear();
  vllhist_2D.clear();
  tthist_2D.clear();
  tthist_alt_2D.clear();
  dbhist_2D.clear();
  gmhist_2D.clear();
  qcdhist_2D.clear();
  ewkwhist_2D.clear();
  ewkzhist_2D.clear();
  dthist_2D.clear();

  if(vnntree) delete vnntree;
  if(vltree) delete vltree;
  if(vlltree) delete vlltree;
  if(tttree) delete tttree;
  if(tttree_alt) delete tttree_alt;
  if(dbtree) delete dbtree;
  if(gmtree) delete gmtree;
  if(qcdtree) delete qcdtree;
  if(ewkwtree) delete ewkwtree;
  if(ewkztree) delete ewkztree;
  if(dttree) delete dttree;
  
  cout << "Templates for the QCD region computed ..." << endl;
}


/////////
void taudatamchist(TFile* outfile,
                   const Category & category,
                   vector<string> observables,
		   vector<string> observables_2D,
                   const double & lumi,
		   const SamplesNLO & nloSamples,
		   const bool & doAlternativeTop    = false,
		   const bool & isHInv              = false,
		   const bool & applyPFWeight       = false) {


  // Files for Znunu,Wlnu, Zll, top, qcd , diboson, signal, data                                                                                                            
  TChain* vnntree = new TChain("tree/tree");
  TChain* vltree  = new TChain("tree/tree");
  TChain* vlltree = new TChain("tree/tree");
  TChain* tttree  = new TChain("tree/tree");
  TChain* tttree_alt = NULL;
  if(doAlternativeTop)
    tttree_alt    = new TChain("tree/tree");
  TChain* dbtree  = new TChain("tree/tree");
  TChain* gmtree  = new TChain("tree/tree");
  TChain* qcdtree = new TChain("tree/tree");
  TChain* ewkwtree = new TChain("tree/tree");
  TChain* ewkztree = new TChain("tree/tree");
  TChain* dttree   = new TChain("tree/tree");

  vnntree->Add((baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/*root").c_str());
  vltree->Add((baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/*root").c_str());
  vlltree->Add((baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/sigfilter/*root").c_str());
  tttree->Add((baseInputTreePath+"/TopAlternative/sigfilter/*root").c_str());
  if(doAlternativeTop and tttree_alt != NULL){
    tttree_alt->Add((baseInputTreePath+"Top/sigfilter/*root").c_str());
  }

  qcdtree->Add((baseInputTreePath+"/QCD/sigfilter/*root").c_str());
  dbtree->Add((baseInputTreePath+"/DiBoson/sigfilter/*root").c_str());
  gmtree->Add((baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/sigfilter/*root").c_str());
  ewkztree->Add((baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/*root").c_str());
  ewkwtree->Add((baseInputTreePath+"/WJetsEWK/sigfilter/*root").c_str());
  dttree->Add((baseInputTreePath+"/MET_reMiniAOD/sigfilter/*root").c_str());

  // make met histograms                                                                                                                                                        
  vector<TH1*> vnnhist;
  vector<TH1*> vlhist;
  vector<TH1*> vllhist;
  vector<TH1*> tthist;
  vector<TH1*> tthist_alt;
  vector<TH1*> dbhist;
  vector<TH1*> gmhist;
  vector<TH1*> qcdhist;
  vector<TH1*> ewkwhist;
  vector<TH1*> ewkzhist;
  vector<TH1*> dthist;

  /////////
  vector<double> bins;
  for(auto obs : observables){

    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable "+obs+" --> please define it"<<endl;

    TH1F* vnnhist_temp = new TH1F(("vnnbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vlhist_temp = new TH1F(("vlbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* vllhist_temp = new TH1F(("vllbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* tthist_temp = new TH1F(("tbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dbhist_temp = new TH1F(("dbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* gmhist_temp = new TH1F(("gbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* qcdhist_temp = new TH1F(("qbkghisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkwhist_temp = new TH1F(("ewkbkgwhisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* ewkzhist_temp = new TH1F(("ewkbkgzhisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
    TH1F* dthist_temp = new TH1F(("datahisttaun_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);

    vnnhist.push_back(dynamic_cast<TH1*>(vnnhist_temp));
    vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
    vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));
    tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
    qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
    dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
    gmhist.push_back(dynamic_cast<TH1*>(gmhist_temp));
    ewkwhist.push_back(dynamic_cast<TH1*>(ewkwhist_temp));
    ewkzhist.push_back(dynamic_cast<TH1*>(ewkzhist_temp));
    dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    
    if(doAlternativeTop){
      TH1F* tthist_alt_temp = new TH1F(("tbkghisttaun_alt_"+obs).c_str(),"",int(bins.size()-1),&bins[0]);
      tthist_alt.push_back(dynamic_cast<TH1*>(tthist_alt_temp));
    }    
  }

  vector<TH2*> vnnhist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vllhist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> tthist_alt_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> gmhist_2D;
  vector<TH2*> qcdhist_2D;
  vector<TH2*> ewkwhist_2D;
  vector<TH2*> ewkzhist_2D;
  vector<TH2*> dthist_2D;


  for(auto obs : observables_2D){

    bin2D bins = selectBinning2D(obs,category);
    if(bins.binX.empty() or bins.binY.empty())
      cout<<"No binning for this observable "+obs+"--> please define it"<<endl;
    
    TH2F* vnnhist_temp = new TH2F(("vnnbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vlhist_temp = new TH2F(("vlbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* vllhist_temp = new TH2F(("vllbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* tthist_temp = new TH2F(("tbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dbhist_temp = new TH2F(("dbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* gmhist_temp = new TH2F(("gbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* qcdhist_temp = new TH2F(("qbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkwhist_temp = new TH2F(("ewkwbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* ewkzhist_temp = new TH2F(("ewkzbkghisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
    TH2F* dthist_temp = new TH2F(("datahisttaun_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

    vnnhist_2D.push_back(dynamic_cast<TH2*>(vnnhist_temp));
    vlhist_2D.push_back(dynamic_cast<TH2*>(vlhist_temp));
    vllhist_2D.push_back(dynamic_cast<TH2*>(vllhist_temp));
    tthist_2D.push_back(dynamic_cast<TH2*>(tthist_temp));
    qcdhist_2D.push_back(dynamic_cast<TH2*>(qcdhist_temp));
    dbhist_2D.push_back(dynamic_cast<TH2*>(dbhist_temp));
    gmhist_2D.push_back(dynamic_cast<TH2*>(gmhist_temp));
    ewkwhist_2D.push_back(dynamic_cast<TH2*>(ewkwhist_temp));
    ewkzhist_2D.push_back(dynamic_cast<TH2*>(ewkzhist_temp));
    dthist_2D.push_back(dynamic_cast<TH2*>(dthist_temp));

    if(doAlternativeTop){
      TH2F* tthist_alt_temp = new TH2F(("tbkghist_alt_"+obs+"_2D").c_str(),"",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);
      tthist_alt_2D.push_back(dynamic_cast<TH2*>(tthist_alt_temp));
    }    
  }

  
  // get k-factors NLO                                                                                                                                                         
  TFile kffile (kfactorFile.c_str());
  TH1*  znlohist = (TH1*) kffile.Get("ZJets_012j_NLO/nominal");
  TH1*  zlohist  = (TH1*) kffile.Get("ZJets_LO/inv_pt");
  TH1*  zewkhist = (TH1*) kffile.Get("EWKcorr/Z");

  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  TH1*  wnlohist = (TH1*) kffile.Get("WJets_012j_NLO/nominal");
  TH1*  wlohist  = (TH1*) kffile.Get("WJets_LO/inv_pt");
  TH1*  wewkhist = (TH1*) kffile.Get("EWKcorr/W");
 
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);

  TH1* anlohist  = (TH1*) kffile.Get("GJets_1j_NLO/nominal_G");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");

  if(aewkhist)
    aewkhist->Divide(anlohist);
  if(anlohist)
    anlohist->Divide(alohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> dyhists;
  vector<TH1*> whists;
  vector<TH1*> ahists;

  // apply EWK and QCD corrections                                                                                                                                              
  float scale = 1;
  if(nloSamples.useZJetsNLO){
    zhists.push_back(zewkhist);
    scale = 3;
  }
  else{
    zhists.push_back(znlohist); 
    zhists.push_back(zewkhist);
  }
  
  if(nloSamples.useDYJetsNLO)
    dyhists.push_back(zewkhist);
  else{
    dyhists.push_back(znlohist); 
    dyhists.push_back(zewkhist);
  }
  
  if(nloSamples.useWJetsNLO)
    whists.push_back(wewkhist);    
  else{
    whists.push_back(wnlohist); 
    whists.push_back(wewkhist);    
  }

  if(nloSamples.usePhotonJetsNLO)
    ahists.push_back(aewkhist);
  else{
    ahists.push_back(anlohist);
    ahists.push_back(aewkhist);
  }


  //////////////////
  bool isWJet = false;
  if(category == Category::monoV)
    isWJet = true;
 
  cout<<"Tau region: Z->nunu sample "<<endl;
  if(not nloSamples.useZJetsNLO)
    makehist4(vnntree,vnnhist,vnnhist_2D,true,Sample::taun,category,false,1.00,lumi,zhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  else
    makehist4(vnntree,vnnhist,vnnhist_2D,true,Sample::taun,category,false,scale,lumi,zhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  cout<<"Tau region: W+jets sample "<<endl;
  makehist4(vltree,vlhist,vlhist_2D,true,Sample::taun,category,false,1.00,lumi,whists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"Tau region: Z+jets sample "<<endl;
  makehist4(vlltree,vllhist,vllhist_2D,true,Sample::taun,category,false,1.00,lumi,dyhists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"Tau region: gamma+jets sample "<<endl;
  makehist4(gmtree,gmhist,gmhist_2D,true,Sample::taun,category,false,1.00,lumi,ahists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"Tau region: ewkw+jets sample "<<endl;
  makehist4(ewkwtree,ewkwhist,ewkwhist_2D,true,Sample::taun,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"Tau region: ewkz+jets sample "<<endl;
  makehist4(ewkztree,ewkzhist,ewkzhist_2D,true,Sample::taun,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight); // temp fix for a wrong xsec
  cout<<"Tau region: TTbar sample "<<endl;
  makehist4(tttree,tthist,tthist_2D,true,Sample::taun,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

    //alternative ttbar             
  if(doAlternativeTop){
    cout<<"Tau region: TTbar alternative sample "<<endl;
    makehist4(tttree_alt,tthist_alt,tthist_alt_2D,true,Sample::taun,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  }

  cout<<"Tau region: Diboson sample "<<endl;
  makehist4(dbtree,dbhist,dbhist_2D,true,Sample::taun,category,isWJet,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  cout<<"Tau region: QCD sample "<<endl;
  makehist4(qcdtree,qcdhist,qcdhist_2D,true,Sample::taun,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);

  // take average of ttbar                                                                                                                                                     
  if(doAlternativeTop){
    for(size_t iHisto = 0; iHisto < tthist.size(); iHisto++)
      makeAverage(tthist.at(iHisto),tthist_alt.at(iHisto));      
    for(size_t iHisto = 0; iHisto < tthist_2D.size(); iHisto++)
      makeAverage(tthist_2D.at(iHisto),tthist_alt_2D.at(iHisto));
  }

  //smooth
  if(doSmoothing){
    for(auto hist : tthist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qcdhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    
    for(auto hist : tthist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : gmhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : dbhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : vllhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : qcdhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkwhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
    for(auto hist : ewkzhist_2D){ if(TString(hist->GetName()).Contains("_met")) smoothEmptyBins(hist,2);}
  }

  // data                                                
  cout<<"Tau region analysis --> loop on data "<<endl;
  makehist4(dttree,dthist,dthist_2D,false,Sample::taun,category,false,1.00,lumi,ehists,"",false,reweightNVTX,0,isHInv,applyPFWeight);
  
  outfile->cd();
  if(not outfile->GetDirectory("Tau"))
    outfile->mkdir("Tau");
  outfile->cd("Tau");
  // store histograms                                                                                                                                                        
  for(auto hist : vnnhist) hist->Write();
  for(auto hist : vlhist) hist->Write();
  for(auto hist : vllhist) hist->Write();
  for(auto hist : tthist) hist->Write();
  for(auto hist : dbhist) hist->Write();
  for(auto hist : gmhist) hist->Write();
  for(auto hist : qcdhist) hist->Write();
  for(auto hist : ewkwhist) hist->Write();
  for(auto hist : ewkzhist) hist->Write();
  for(auto hist : dthist) hist->Write();

  // store hist_2Dograms                                                                                                                                                     
  outfile->cd();
  outfile->cd("Tau");
  for(auto hist_2D : vnnhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vlhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : vllhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : gmhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : tthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dbhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : qcdhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkwhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : ewkzhist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }
  for(auto hist_2D : dthist_2D){ TH1* temp = unroll2DHistograms(hist_2D); temp->Write(); }

  outfile->cd();

  for(auto hist : vnnhist) delete hist;
  for(auto hist : vlhist) delete hist;
  for(auto hist : vllhist) delete hist;
  for(auto hist : tthist) delete hist;
  for(auto hist : dbhist) delete hist;
  for(auto hist : gmhist) delete hist;
  for(auto hist : qcdhist) delete hist;
  for(auto hist : ewkwhist) delete hist;
  for(auto hist : ewkzhist) delete hist;
  for(auto hist : dthist) delete hist;

  vnnhist.clear();
  vlhist.clear();
  vllhist.clear();
  tthist.clear();
  dbhist.clear();
  gmhist.clear();
  qcdhist.clear();
  ewkwhist.clear();
  ewkzhist.clear();
  dthist.clear();

  for(auto hist : vnnhist_2D) delete hist;
  for(auto hist : vlhist_2D) delete hist;
  for(auto hist : vllhist_2D) delete hist;
  for(auto hist : tthist_2D) delete hist;
  for(auto hist : dbhist_2D) delete hist;
  for(auto hist : gmhist_2D) delete hist;
  for(auto hist : qcdhist_2D) delete hist;
  for(auto hist : ewkwhist_2D) delete hist;
  for(auto hist : ewkzhist_2D) delete hist;
  for(auto hist : dthist_2D) delete hist;

  vnnhist_2D.clear();
  vlhist_2D.clear();
  vllhist_2D.clear();
  tthist_2D.clear();
  tthist_alt_2D.clear();
  dbhist_2D.clear();
  gmhist_2D.clear();
  qcdhist_2D.clear();
  ewkwhist_2D.clear();
  ewkzhist_2D.clear();
  dthist_2D.clear();

  if(vnntree) delete vnntree;
  if(vltree) delete vltree;
  if(vlltree) delete vlltree;
  if(tttree) delete tttree;
  if(tttree_alt) delete tttree_alt;
  if(dbtree) delete dbtree;
  if(gmtree) delete gmtree;
  if(qcdtree) delete qcdtree;
  if(ewkwtree) delete ewkwtree;
  if(ewkztree) delete ewkztree;
  if(dttree) delete dttree;
  
  cout << "Templates for the Tau region computed ..." << endl;
}
