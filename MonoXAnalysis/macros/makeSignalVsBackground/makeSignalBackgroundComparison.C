//monoV_signalAnalysis("/store/ciMassaf/user/rgerosa/MONOJET_ANALYSIS/Production-14-1-2016/",1,"Pseudoscalar","Plots_MonoW_Pseudoscalar")
#include <iostream>
#include <vector>
#include <string>
#include <string>
#include <fstream>
#include <unordered_map>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"


vector<double> bins_monoJ_met         = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
vector<double> bins_monoV_met         = {250.,300.,350.,400.,500.,600.,750.,1000.};
vector<double> bins_monoJ_dphiJJ      = {-0.1,0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.14};
//vector<double> bins_monoJ_dphiJJ      = {-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.14};
vector<double> bins_monoV_dphiJJ      = {-0.1,0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,3.14};
vector<double> bins_monoJ_dRJJ        = {-0.1,0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6,6.5,7.,7.5};
vector<double> bins_monoV_dRJJ        = {-0.1,0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6,6.5,7.,7.5};

vector<double> bins_monoJ_mT          = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1350,1550,1750,2000};
vector<double> bins_monoV_mT          = {200.,250.,300.,350.,400.,500.,600.,1000,1250,1500,1750,2000};
vector<double> bins_monoJ_HT          = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
vector<double> bins_monoV_HT          = {200.,300.,400.,500,600.,700.,850,950,1050,1250,1650,2000};
vector<double> bins_monoV_njet        = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
vector<double> bins_monoJ_njet        = {0.,1.,2.,3.,4.,5.,6.,7.,8.};

vector<double> bins_monoV_mpr        = {65.,67.5,70.,72.5,75.,77.5,80.,82.5,85.,87.5,90.,92.5,95.,97.5,100.,102.5,105.};
vector<double> bins_monoJ_mpr        = {0.,4.,8.,12.,16.,20.,24.,28.,32.,36.,40.,44.,48.,52.,56.,60.,64.,68.,72.,76.,80.,84.,88.,92.,96.,100.,104.,110.,116.,125.,135.};
vector<double> bins_monoV_QGL        = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};
vector<double> bins_monoJ_QGL        = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};
vector<double> bins_monoJ_jetPt      = {100.,120.,140.,160.,180.,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
vector<double> bins_monoV_jetPt      = {250,280,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
vector<double> bins_monoV_tau2tau1   = {0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
vector<double> bins_monoJ_tau2tau1   = {0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};

vector<int> color = {2,4,6,7,12,3,9,46};

using namespace std;
void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList_MV, vector<TH1F*> histoList_MJ,
		    TH1F* histBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTile);
void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList_MV, 
		    TH1F* histBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTile);

// make signal efficiency studies
void makeSignalBackgroundComparison(string baseInputPath,    // base path with all the directories on the local pc
				    string outputDirectory,  // output directory for plots
				    bool   isHiggsInvisible,
				    string interactionModel, // interaction 				    
				    string bosonMode, // MonoW or MonoZ depending on which mono-V signal to display
				    int    category, // 3 no monoV cuts, 1 mono-jet channel, 2 mono-V channel
				    vector<string> Mphi,
				    vector<string> Mchi,
				    bool   displayOnlyMonoV
				    ){

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);
  gStyle->SetOptTitle(0);

  //create output directories
  system(("mkdir -p "+outputDirectory).c_str());

  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1*   puhist = (TH1*) pufile->Get("puhist");

  // trigger efficiency for met trigger                                                                                                                                         
  TFile* trmfile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/mettrigSF.root");
  TH1*   trmhist = (TH1*) trmfile->Get("mettrigSF");


  //////
  TCanvas *cCanvas = new TCanvas("cCanvas","",180,52,550,550);
  cCanvas->SetTicks();
  cCanvas->SetFillColor(0);
  cCanvas->SetBorderMode(0);
  cCanvas->SetBorderSize(2);
  cCanvas->SetRightMargin(0.05);
  cCanvas->SetBottomMargin(0.12);
  cCanvas->SetFrameBorderMode(0);
  
  TLegend* legend = new TLegend(0.18,0.68,0.7,0.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.031);
  legend->SetTextFont(42);

  float lumi = 2.30;
  
  // take signal region files for signal
  vector<TFile*> fileMonoV_Sig;  
  vector<TFile*> fileMonoJet_Sig;  
  vector<pair<string,TTree*> > treeMonoV_Sig;  
  vector<pair<string,TTree*> > treeMonoJet_Sig;  

  string jetMode;
  ifstream infile;
  string line;

  if(not isHiggsInvisible){

    if(interactionModel == "Vector" or interactionModel == "Axial")
      jetMode = "DMV";
    else if(interactionModel == "Scalar" or interactionModel == "PseudoScalar")
      jetMode = "DMS";

    system(("ls "+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/ > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);
	if(line == "") continue;	
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	    cout<<baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/"+line<<endl;
	    fileMonoV_Sig.push_back(TFile::Open((baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/"+line).c_str()));
	    if(fileMonoV_Sig.back() and not fileMonoV_Sig.back()->IsZombie())
	    treeMonoV_Sig.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_Sig.back()->Get("tree/tree")));
	  }
	}
      }
    }

    infile.close();

    system(("ls "+baseInputPath+"/"+jetMode+"_"+interactionModel+"/sigfilter/ > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);
	if(line == "") continue;	
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	    cout<<baseInputPath+"/"+jetMode+"_"+interactionModel+"/sigfilter/"+line<<endl;
	    fileMonoJet_Sig.push_back(TFile::Open((baseInputPath+"/"+jetMode+"_"+interactionModel+"/sigfilter/"+line).c_str()));
	    if(fileMonoJet_Sig.back() and not fileMonoJet_Sig.back()->IsZombie())
	    treeMonoJet_Sig.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoJet_Sig.back()->Get("tree/tree")));
	  }
	}
      }
    }
    infile.close();
  }

  else{
    
    if (bosonMode == "MonoW")
      bosonMode = "WH";
    else if (bosonMode == "MonoZ")
      bosonMode = "ZH";
    jetMode = "ggH";
    
    system(("ls "+baseInputPath+"/HiggsInvisible/sigfilter/ | grep "+bosonMode+" > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
        getline(infile,line);
        if(line == "") continue;
        for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
          if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
            cout<<"found monoV file "<<baseInputPath+"/HiggsInvisible/sigfilter/"+line<<endl;
            fileMonoV_Sig.push_back(TFile::Open((baseInputPath+"/HiggsInvisible/sigfilter/"+line).c_str()));
            if(fileMonoV_Sig.back() and not fileMonoV_Sig.back()->IsZombie())
	      treeMonoV_Sig.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_Sig.back()->Get("tree/tree")));
          }
        }
      }
    }

    infile.close();

    system(("ls "+baseInputPath+"/HiggsInvisible/sigfilter/ | grep GluGlu_HToInvisible > file_list.txt").c_str());    
    infile.open("file_list.txt",ifstream::in);

    if(infile.is_open()){
      while(!infile.eof()){
        getline(infile,line);
        if(line == "") continue;
        for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
          if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
            cout<<"found ggH file : "+baseInputPath+"/HiggsInvisible//sigfilter/"+line<<endl;
            fileMonoJet_Sig.push_back(TFile::Open((baseInputPath+"/HiggsInvisible/sigfilter/"+line).c_str()));
            if(fileMonoJet_Sig.back() and not fileMonoJet_Sig.back()->IsZombie())
	      treeMonoJet_Sig.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoJet_Sig.back()->Get("tree/tree")));
          }
        }
      }
    }
    infile.close();    
  }


  //Background Znunu
  TChain* backgroundZnunu = new TChain("tree/tree");
  TFile* file1 = TFile::Open((baseInputPath+"/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root").c_str());
  backgroundZnunu->Add(file1->GetName());

  string cut;
  if (category == 2)
    cut = "(hltmet90 || hltmet120 || hltmetwithmu120 || hltmetwithmu170 || hltmetwithmu300 || hltmetwithmu90) && njets >= 1 && combinejetpt[0] > 100 && abs(combinejeteta[0]) < 2.5 && nbjetslowpt == 0 && combinejetCHfrac[0] > 0.1 && combinejetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200 && nmuons == 0 && nphotons == 0 && ntaus == 0 && nelectrons == 0 && fabs(boostedJeteta[0]) < 2.4 &&  boostedJetpt[0] > 250 && prunedJetm[0] > 65 && prunedJetm[0] < 105 && boostedJettau2[0]/boostedJettau1[0] < 0.6 && t1pfmet > 250";
  else if(category == 1)
    cut = "(hltmet90 || hltmet120 || hltmetwithmu120 || hltmetwithmu170 || hltmetwithmu300 || hltmetwithmu90) && njets >= 1 && combinejetpt[0] > 100 && abs(combinejeteta[0]) < 2.5 && nbjetslowpt == 0 && combinejetCHfrac[0] > 0.1 && combinejetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200 && nmuons == 0 && nphotons == 0 && ntaus == 0 && nelectrons == 0 && (fabs(boostedJeteta[0]) > 2.4 || boostedJetpt[0] < 250 || prunedJetm[0] < 65 || prunedJetm[0] > 105 || boostedJettau2[0]/boostedJettau1[0] > 0.6 || t1pfmet < 250)";
  else if(category == 0)
    cut = "(hltmet90 || hltmet120 || hltmetwithmu120 || hltmetwithmu170 || hltmetwithmu300 || hltmetwithmu90) && njets >= 1 && combinejetpt[0] > 100 && abs(combinejeteta[0]) < 2.5 && nbjetslowpt == 0 && combinejetCHfrac[0] > 0.1 && combinejetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200 && nmuons == 0 && nphotons == 0 && ntaus == 0 && nelectrons == 0 ";
  
  /// now look at substructure in the signal region only
  vector<TH1F*> boostedJetPt_monoV;
  vector<TH1F*> boostedJetEta_monoV;
  vector<TH1F*> boostedJetTau2Tau1_monoV;
  vector<TH1F*> boostedPrunedJetm_monoV;
  vector<TH1F*> leadingJetPt_monoV;  
  vector<TH1F*> leadingJetEta_monoV;
  vector<TH1F*> njets_monoV;
  vector<TH1F*> njets_fwd_monoV;
  vector<TH1F*> ht_monoV;
  vector<TH1F*> mT_monoV;
  vector<TH1F*> met_monoV;
  vector<TH1F*> dphijj_monoV;
  vector<TH1F*> mindphijj_monoV;
  vector<TH1F*> mindphij1j_monoV;
  vector<TH1F*> dRjj_monoV;
  vector<TH1F*> mindRj1j_monoV;
  vector<TH1F*> leadingJetQGL_monoV;

  vector<TH1F*> boostedJetPt_monoJet;
  vector<TH1F*> boostedJetEta_monoJet;
  vector<TH1F*> boostedJetTau2Tau1_monoJet;
  vector<TH1F*> boostedPrunedJetm_monoJet;
  vector<TH1F*> leadingJetPt_monoJet;  
  vector<TH1F*> leadingJetEta_monoJet;
  vector<TH1F*> njets_monoJet;
  vector<TH1F*> njets_fwd_monoJet;
  vector<TH1F*> ht_monoJet;
  vector<TH1F*> mT_monoJet;
  vector<TH1F*> met_monoJet;
  vector<TH1F*> dphijj_monoJet;
  vector<TH1F*> mindphijj_monoJet;
  vector<TH1F*> mindphij1j_monoJet;
  vector<TH1F*> dRjj_monoJet;
  vector<TH1F*> mindRj1j_monoJet;
  vector<TH1F*> leadingJetQGL_monoJet;

  cout<<"Run signal analysis "<<endl;

  for(size_t iMass = 0; iMass < treeMonoV_Sig.size(); iMass++){

    string name = treeMonoV_Sig.at(iMass).first;
    TTreeFormula* formula_MV = new TTreeFormula((bosonMode+"_"+interactionModel+"_"+name+"_formula").c_str(),
						cut.c_str(),treeMonoV_Sig.at(iMass).second);
    // allocate all the histograms
    cout<<"histogram name "<<bosonMode+"_"+interactionModel+"_"+name<<endl;
    leadingJetEta_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetEta").c_str(),(bosonMode+"_"+interactionModel+"_"+name).c_str(),20,-2.5,2.5));  
    boostedJetEta_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedJetEta").c_str(),(bosonMode+"_"+interactionModel+"_"+name).c_str(),20,-2.5,2.5));  

    if(category == 2){
      boostedJetPt_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedJetPt").c_str(),
					    (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_jetPt.size()-1,&bins_monoV_jetPt[0]));
      boostedJetTau2Tau1_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedJetTau2Tau1").c_str(),
						  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_tau2tau1.size()-1,&bins_monoV_tau2tau1[0]));    

      boostedPrunedJetm_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedPrunedJetm").c_str(),
						 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_mpr.size()-1,&bins_monoV_mpr[0]));
      leadingJetPt_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetPt").c_str(),
					    (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_jetPt.size()-1,&bins_monoV_jetPt[0]));
      njets_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_njets").c_str(),
				     (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_njet.size()-1,&bins_monoV_njet[0]));
      njets_fwd_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_njets_fwd").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_njet.size()-1,&bins_monoV_njet[0]));
      ht_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_ht").c_str(),
				  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_HT.size()-1,&bins_monoV_HT[0]));
      mT_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mT").c_str(),				
				  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_mT.size()-1,&bins_monoV_mT[0]));
      met_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_met").c_str(),
				   (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_met.size()-1,&bins_monoV_met[0]));
      dphijj_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_dphijj").c_str(),
				      (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dphiJJ.size()-1,&bins_monoV_dphiJJ[0]));
      mindphijj_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mindphijj").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dphiJJ.size()-1,&bins_monoV_dphiJJ[0]));
      mindphij1j_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mindphij1j").c_str(),
					  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dphiJJ.size()-1,&bins_monoV_dphiJJ[0]));
      leadingJetQGL_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetQGL").c_str(),
					     (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_QGL.size()-1,&bins_monoV_QGL[0]));

      dRjj_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_dRjj").c_str(),
				    (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dRJJ.size()-1,&bins_monoV_dRJJ[0]));
      mindRj1j_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mindRj1j").c_str(),
					(bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dRJJ.size()-1,&bins_monoV_dRJJ[0]));
    }
    else{
      
      boostedJetPt_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedJetPt").c_str(),
					    (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]));      
      boostedJetTau2Tau1_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedJetTau2Tau1").c_str(),
						  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_tau2tau1.size()-1,&bins_monoJ_tau2tau1[0]));    
      boostedPrunedJetm_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedPrunedJetm").c_str(),
						 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_mpr.size()-1,&bins_monoJ_mpr[0]));
      leadingJetPt_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetPt").c_str(),
					    (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]));
      njets_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_njets").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]));    
      njets_fwd_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_njets_fwd").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]));    
      ht_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_ht").c_str(),
				  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_HT.size()-1,&bins_monoJ_HT[0]));
      mT_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mT").c_str(),				
				  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_mT.size()-1,&bins_monoJ_mT[0]));
      met_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_met").c_str(),
				   (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_met.size()-1,&bins_monoJ_met[0]));
      dphijj_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_dphijj").c_str(),
				      (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      mindphijj_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mindphijj").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      mindphij1j_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mindphij1j").c_str(),
					  (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));

      leadingJetQGL_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetQGL").c_str(),
					     (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_QGL.size()-1,&bins_monoJ_QGL[0]));

      dRjj_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_dRjj").c_str(),
				    (bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dRJJ.size()-1,&bins_monoJ_dRJJ[0]));
      mindRj1j_monoV.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_mindRj1j").c_str(),
					(bosonMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dRJJ.size()-1,&bins_monoJ_dRJJ[0]));
    }
  
    boostedJetPt_monoV.back()->Sumw2();
    boostedJetEta_monoV.back()->Sumw2();
    boostedJetTau2Tau1_monoV.back()->Sumw2();
    boostedPrunedJetm_monoV.back()->Sumw2();
    leadingJetPt_monoV.back()->Sumw2();
    leadingJetEta_monoV.back()->Sumw2();
    njets_monoV.back()->Sumw2();
    njets_fwd_monoV.back()->Sumw2();
    ht_monoV.back()->Sumw2();
    mT_monoV.back()->Sumw2();
    met_monoV.back()->Sumw2();
    dphijj_monoV.back()->Sumw2();
    mindphijj_monoV.back()->Sumw2();
    mindphij1j_monoV.back()->Sumw2();
    leadingJetQGL_monoV.back()->Sumw2();
    dRjj_monoV.back()->Sumw2();
    mindRj1j_monoV.back()->Sumw2();
    
    TTreeReader myReader_monoV(treeMonoV_Sig.at(iMass).second);
    TTreeReaderValue<vector<double> > boostedJetpt_MV    (myReader_monoV,"boostedJetpt");
    TTreeReaderValue<vector<double> > boostedJeteta_MV   (myReader_monoV,"boostedJeteta");
    TTreeReaderValue<vector<double> > boostedJetphi_MV   (myReader_monoV,"boostedJetphi");
    TTreeReaderValue<vector<double> > boostedJetm_MV     (myReader_monoV,"boostedJetm");
    TTreeReaderValue<vector<double> > boostedJettau2_MV  (myReader_monoV,"boostedJettau2");
    TTreeReaderValue<vector<double> > boostedJettau1_MV  (myReader_monoV,"boostedJettau1");
    TTreeReaderValue<vector<double> > boostedJetm_pr_MV  (myReader_monoV,"prunedJetm");

    TTreeReaderValue<vector<double> > jetpt_MV    (myReader_monoV,"combinejetpt");
    TTreeReaderValue<vector<double> > jeteta_MV   (myReader_monoV,"combinejeteta");
    TTreeReaderValue<vector<double> > jetphi_MV   (myReader_monoV,"combinejetphi");
    TTreeReaderValue<vector<double> > jetQGL_MV   (myReader_monoV,"combinejetQGL");

    TTreeReaderValue<unsigned int> njets_MV  (myReader_monoV,"njets");
    TTreeReaderValue<vector<double> > fwdjetpt_MV  (myReader_monoV,"forwardjetpt");
    TTreeReaderValue<double> ht_MV  (myReader_monoV,"ht");
    TTreeReaderValue<double> met_MV  (myReader_monoV,"t1pfmet");
    TTreeReaderValue<double> metphi_MV  (myReader_monoV,"t1pfmetphi");

    //    TTreeReaderValue<double> xsec_MV         (myReader_monoV,"xsec");
    //    TTreeReaderValue<double> wgt_MV          (myReader_monoV,"wgt");
    //    TTreeReaderValue<double> wgtbtag_MV      (myReader_monoV,"wgtbtag");
    //    TTreeReaderValue<double> wgtsum_MV       (myReader_monoV,"wgtsum");
    TTreeReaderValue<unsigned int> nvtx_MV   (myReader_monoV,"nvtx");

    int iEvent = 0;
    cout<<"Run on monoV events "<<treeMonoV_Sig.at(iMass).second->GetEntries()<<endl;

    while(myReader_monoV.Next()){
      treeMonoV_Sig.at(iMass).second->GetEntry(iEvent); 
      cout.flush();
      if(iEvent %5000 == 0) cout<<"\r"<<"iEvent "<<100*float(iEvent)/treeMonoV_Sig.at(iMass).second->GetEntries()<<" %";
      iEvent++;
      if(not formula_MV->EvalInstance()) continue;

      float weight = (puhist->GetBinContent(*nvtx_MV))*(trmhist->GetBinContent(trmhist->FindBin(*met_MV)));

      if(boostedJetpt_MV->size()>0)
	boostedJetPt_monoV.back()->Fill(boostedJetpt_MV->at(0),weight);
      if(boostedJeteta_MV->size()>0)
	boostedJetEta_monoV.back()->Fill(boostedJeteta_MV->at(0),weight);
      if(boostedJettau2_MV->size()>0 and boostedJettau1_MV->size()>0 )
	boostedJetTau2Tau1_monoV.back()->Fill(boostedJettau2_MV->at(0)/boostedJettau1_MV->at(0),weight);
      if(boostedJetm_pr_MV->size()>0)
	boostedPrunedJetm_monoV.back()->Fill(boostedJetm_pr_MV->at(0),weight);
      if(jetpt_MV->size()>0)
	leadingJetPt_monoV.back()->Fill(jetpt_MV->at(0),weight);
      if(jeteta_MV->size()>0)
	leadingJetEta_monoV.back()->Fill(jeteta_MV->at(0),weight);
      if(jetQGL_MV->size()>0)
	leadingJetQGL_monoV.back()->Fill(jetQGL_MV->at(0),weight);

      njets_monoV.back()->Fill(*njets_MV,weight);
      njets_fwd_monoV.back()->Fill(*njets_MV+fwdjetpt_MV->size(),weight);
      ht_monoV.back()->Fill(*ht_MV,weight);
      met_monoV.back()->Fill(*met_MV,weight);

      if(jetpt_MV->size()>0){
	float deltaPhi = fabs(jetphi_MV->at(0)-*metphi_MV);
	if(deltaPhi > TMath::Pi())
	  deltaPhi = fabs(2*TMath::Pi() - deltaPhi);
	mT_monoV.back()->Fill(sqrt(2*jetpt_MV->at(0)*(*met_MV)*(1-cos(deltaPhi))),weight);
      }


      if(jetphi_MV->size() < 2){
	dphijj_monoV.back()->Fill(dphijj_monoV.back()->GetBinCenter(1),weight);
	dRjj_monoV.back()->Fill(dRjj_monoV.back()->GetBinCenter(1),weight);
      }
      else{
	if(fabs(jetphi_MV->at(0)-jetphi_MV->at(1)) > TMath::Pi()){
	  dphijj_monoV.back()->Fill(2*TMath::Pi()-fabs(jetphi_MV->at(0)-jetphi_MV->at(1)),weight);
	  dRjj_monoV.back()->Fill(sqrt(TMath::Power(2*TMath::Pi()-fabs(jetphi_MV->at(0)-jetphi_MV->at(1)),2)+TMath::Power(jeteta_MV->at(0)-jeteta_MV->at(1),2)),weight);
	}
	else{
	  dphijj_monoV.back()->Fill(fabs(jetphi_MV->at(0)-jetphi_MV->at(1)),weight);
	  dRjj_monoV.back()->Fill(sqrt(TMath::Power(fabs(jetphi_MV->at(0)-jetphi_MV->at(1)),2)+TMath::Power(jeteta_MV->at(0)-jeteta_MV->at(1),2)),weight);
	}
	
      }

      if(jetphi_MV->size() < 2){
	mindphijj_monoV.back()->Fill(mindphijj_monoV.back()->GetBinCenter(1),weight);
	mindphij1j_monoV.back()->Fill(mindphij1j_monoV.back()->GetBinCenter(1),weight);
	mindRj1j_monoV.back()->Fill(dRjj_monoV.back()->GetBinCenter(1),weight);
      }
      else{
	float minDphijj  = 2*TMath::Pi();
	float minDphij1j = 2*TMath::Pi();
	bool  isfound_jj  = false;
	bool  isfound_j1j = false;
	for(size_t ijet = 0 ; ijet < jetphi_MV->size(); ijet++){
	  for(size_t jjet = ijet+1 ; jjet < jetphi_MV->size(); jjet++){
	    if(jetpt_MV->at(ijet) < 30 or jetpt_MV->at(jjet) < 30) continue; 
	    float deltaPhi = fabs(jetphi_MV->at(ijet)-jetphi_MV->at(jjet));
	    if(deltaPhi > TMath::Pi())
	      deltaPhi = 2*TMath::Pi() - deltaPhi;
	    if(deltaPhi > 0 and deltaPhi < minDphijj){
	      minDphijj  = deltaPhi;
	      isfound_jj = true;
	    }
	  }
	}
	
	for(size_t jjet = 1 ; jjet < jetphi_MV->size(); jjet++){
	  float deltaPhi = fabs(jetphi_MV->at(0)-jetphi_MV->at(jjet));
	  if(jetpt_MV->at(0) < 30 or jetpt_MV->at(jjet) < 30) continue;
	  if(deltaPhi > TMath::Pi())
	    deltaPhi = 2*TMath::Pi() - deltaPhi;
	  if(deltaPhi > 0 and deltaPhi < minDphij1j){
	    minDphij1j  = deltaPhi;
	    isfound_j1j = true;
	  }
	}
      
      
	if(isfound_jj)
	  mindphijj_monoV.back()->Fill(minDphijj,weight);
	else 
	  mindphijj_monoV.back()->Fill(mindphijj_monoV.back()->GetBinCenter(1),weight);

	if(isfound_j1j){
	  mindphij1j_monoV.back()->Fill(minDphij1j,weight);
	  mindRj1j_monoV.back()->Fill(sqrt(minDphij1j*minDphij1j+TMath::Power(jeteta_MV->at(0)-jeteta_MV->at(1),2)),weight);
	}
	else{
	  mindphij1j_monoV.back()->Fill(mindphij1j_monoV.back()->GetBinCenter(1),weight);
	  mindRj1j_monoV.back()->Fill(mindphij1j_monoV.back()->GetBinCenter(1),weight);
	}
      }
    }
    std::cout<<std::endl;
    ////
    leadingJetEta_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_leadingJetEta").c_str(),(jetMode+"_"+interactionModel+"_"+name).c_str(),20,-2.5,2.5));  
    boostedJetEta_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedJetEta").c_str(),(jetMode+"_"+interactionModel+"_"+name).c_str(),20,-2.5,2.5));  
    
    if(category == 2){
      boostedJetPt_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedJetPt").c_str(),
					    (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]));
      boostedJetTau2Tau1_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedJetTau2Tau1").c_str(),
						  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_tau2tau1.size()-1,&bins_monoJ_tau2tau1[0]));    

      boostedPrunedJetm_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedPrunedJetm").c_str(),
						 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_mpr.size()-1,&bins_monoJ_mpr[0]));
      leadingJetPt_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_leadingJetPt").c_str(),
					    (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]));
      njets_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_njets").c_str(),
				     (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]));
      njets_fwd_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_njets_fwd").c_str(),
					 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]));
      ht_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_ht").c_str(),
				  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_HT.size()-1,&bins_monoJ_HT[0]));
      mT_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mT").c_str(),				
				  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_mT.size()-1,&bins_monoJ_mT[0]));
      met_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_met").c_str(),
				   (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_met.size()-1,&bins_monoJ_met[0]));
      dphijj_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_dphijj").c_str(),
				      (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      mindphijj_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mindphijj").c_str(),
					 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      mindphij1j_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mindphij1j").c_str(),
					  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      leadingJetQGL_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_leadingJetQGL").c_str(),
					     (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_QGL.size()-1,&bins_monoJ_QGL[0]));

      dRjj_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_dRjj").c_str(),
                                    (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dRJJ.size()-1,&bins_monoV_dRJJ[0]));
      mindRj1j_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mindRj1j").c_str(),
                                        (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoV_dRJJ.size()-1,&bins_monoV_dRJJ[0]));

    }
 
    else{
      
      boostedJetPt_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedJetPt").c_str(),
					      (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]));      
      boostedJetTau2Tau1_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedJetTau2Tau1").c_str(),
						  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_tau2tau1.size()-1,&bins_monoJ_tau2tau1[0]));    
      boostedPrunedJetm_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_boostedPrunedJetm").c_str(),
						 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_mpr.size()-1,&bins_monoJ_mpr[0]));
      leadingJetPt_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_leadingJetPt").c_str(),
					    (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]));
      njets_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_njets").c_str(),
					 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]));    
      njets_fwd_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_njets_fwd").c_str(),
					 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]));    
      ht_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_ht").c_str(),
				  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_HT.size()-1,&bins_monoJ_HT[0]));
      mT_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mT").c_str(),				
				  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_mT.size()-1,&bins_monoJ_mT[0]));
      met_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_met").c_str(),
				   (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_met.size()-1,&bins_monoJ_met[0]));
      dphijj_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_dphijj").c_str(),
				      (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      mindphijj_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mindphijj").c_str(),
					 (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));
      mindphij1j_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mindphij1j").c_str(),
					  (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]));

      leadingJetQGL_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_leadingJetQGL").c_str(),
					     (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_QGL.size()-1,&bins_monoJ_QGL[0]));
      
      dRjj_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_dRjj").c_str(),
                                    (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dRJJ.size()-1,&bins_monoJ_dRJJ[0]));
      mindRj1j_monoJet.push_back(new TH1F((jetMode+"_"+interactionModel+"_"+name+"_mindRj1j").c_str(),
                                        (jetMode+"_"+interactionModel+"_"+name).c_str(),bins_monoJ_dRJJ.size()-1,&bins_monoJ_dRJJ[0]));
    }
     
    
    boostedJetPt_monoJet.back()->Sumw2();
    boostedJetEta_monoJet.back()->Sumw2();
    boostedJetTau2Tau1_monoJet.back()->Sumw2();
    boostedPrunedJetm_monoJet.back()->Sumw2();
    leadingJetPt_monoJet.back()->Sumw2();
    leadingJetEta_monoJet.back()->Sumw2();
    njets_monoJet.back()->Sumw2();
    njets_fwd_monoJet.back()->Sumw2();
    ht_monoJet.back()->Sumw2();
    mT_monoJet.back()->Sumw2();
    met_monoJet.back()->Sumw2();
    dphijj_monoJet.back()->Sumw2();
    mindphijj_monoJet.back()->Sumw2();
    mindphij1j_monoJet.back()->Sumw2();
    leadingJetQGL_monoJet.back()->Sumw2();
    dRjj_monoJet.back()->Sumw2();
    mindRj1j_monoJet.back()->Sumw2();

    if(not displayOnlyMonoV){
      
      TTreeFormula* formula_MJ = new TTreeFormula((jetMode+"_"+interactionModel+"_"+name+"_formula").c_str(),
						  cut.c_str(),treeMonoJet_Sig.at(iMass).second);
    
      TTreeReader myReader_monoJet(treeMonoJet_Sig.at(iMass).second);
      TTreeReaderValue<vector<double> > boostedJetpt_MJ    (myReader_monoJet,"boostedJetpt");
      TTreeReaderValue<vector<double> > boostedJeteta_MJ   (myReader_monoJet,"boostedJeteta");
      TTreeReaderValue<vector<double> > boostedJetphi_MJ   (myReader_monoJet,"boostedJetphi");
      TTreeReaderValue<vector<double> > boostedJetm_MJ     (myReader_monoJet,"boostedJetm");
      TTreeReaderValue<vector<double> > boostedJettau2_MJ  (myReader_monoJet,"boostedJettau2");
      TTreeReaderValue<vector<double> > boostedJettau1_MJ  (myReader_monoJet,"boostedJettau1");
      TTreeReaderValue<vector<double> > boostedJetm_pr_MJ  (myReader_monoJet,"prunedJetm");
      TTreeReaderValue<vector<double> > jetpt_MJ    (myReader_monoJet,"combinejetpt");
      TTreeReaderValue<vector<double> > jeteta_MJ   (myReader_monoJet,"combinejeteta");
      TTreeReaderValue<vector<double> > jetphi_MJ   (myReader_monoJet,"combinejetphi");
      TTreeReaderValue<vector<double> > jetQGL_MJ   (myReader_monoJet,"combinejetQGL");
      TTreeReaderValue<unsigned int> njets_MJ  (myReader_monoJet,"njets");
      TTreeReaderValue<vector<double> > fwdjetpt_MJ  (myReader_monoJet,"forwardjetpt");
      TTreeReaderValue<double> ht_MJ  (myReader_monoJet,"ht");
      TTreeReaderValue<double> met_MJ  (myReader_monoJet,"t1pfmet");
      TTreeReaderValue<double> metphi_MJ  (myReader_monoJet,"t1pfmetphi");
      //      TTreeReaderValue<double> xsec_MJ         (myReader_monoJet,"xsec");
      //      TTreeReaderValue<double> wgt_MJ          (myReader_monoJet,"wgt");
      //      TTreeReaderValue<double> wgtbtag_MJ      (myReader_monoJet,"wgtbtag");
      //      TTreeReaderValue<double> wgtsum_MJ       (myReader_monoJet,"wgtsum");
      TTreeReaderValue<unsigned int> nvtx_MJ   (myReader_monoJet,"nvtx");
      
      iEvent = 0;
      cout<<"Run on mono-Jet events "<<treeMonoJet_Sig.at(iMass).second->GetEntries()<<endl;
      while(myReader_monoJet.Next()){
	
	float weight = (puhist->GetBinContent(*nvtx_MJ))*(trmhist->GetBinContent(trmhist->FindBin(*met_MJ)));
	
	treeMonoJet_Sig.at(iMass).second->GetEntry(iEvent);
	cout.flush();
	if(iEvent %5000 == 0) cout<<"\r"<<"iEvent "<<100*float(iEvent)/treeMonoJet_Sig.at(iMass).second->GetEntries()<<" %";
	iEvent++;
	if(not formula_MJ->EvalInstance()) continue;
	
	if(boostedJetpt_MJ->size()>0)
	  boostedJetPt_monoJet.back()->Fill(boostedJetpt_MJ->at(0),weight);
	if(boostedJeteta_MJ->size()>0)
	  boostedJetEta_monoJet.back()->Fill(boostedJeteta_MJ->at(0),weight);
	if(boostedJettau2_MJ->size()>0 and boostedJettau1_MJ->size()>0 )
	  boostedJetTau2Tau1_monoJet.back()->Fill(boostedJettau2_MJ->at(0)/boostedJettau1_MJ->at(0),weight);
	if(boostedJetm_pr_MJ->size()>0)
	  boostedPrunedJetm_monoJet.back()->Fill(boostedJetm_pr_MJ->at(0),weight);
	if(jetpt_MJ->size()>0)
	  leadingJetPt_monoJet.back()->Fill(jetpt_MJ->at(0),weight);
	if(jeteta_MJ->size()>0)
	  leadingJetEta_monoJet.back()->Fill(jeteta_MJ->at(0),weight);
	if(jetQGL_MJ->size()>0)
	  leadingJetQGL_monoJet.back()->Fill(jetQGL_MJ->at(0),weight);
	
	njets_monoJet.back()->Fill(*njets_MJ,weight);
	njets_fwd_monoJet.back()->Fill(*njets_MJ+fwdjetpt_MJ->size(),weight);
	ht_monoJet.back()->Fill(*ht_MJ,weight);
	met_monoJet.back()->Fill(*met_MJ,weight);
	
	if(jetpt_MJ->size()>0){
	  float deltaPhi = fabs(jetphi_MJ->at(0)-*metphi_MJ);
	  if(deltaPhi > TMath::Pi())
	    deltaPhi = fabs(2*TMath::Pi() - deltaPhi);
	  mT_monoJet.back()->Fill(sqrt(2*jetpt_MJ->at(0)*(*met_MJ)*(1-cos(deltaPhi))),weight);
	}
	//////
	if(jetphi_MJ->size() < 2){
	  dphijj_monoJet.back()->Fill(dphijj_monoJet.back()->GetBinCenter(1),weight);
	  dRjj_monoJet.back()->Fill(dRjj_monoJet.back()->GetBinCenter(1),weight);
	}
	else{
	  if(fabs(jetphi_MJ->at(0)-jetphi_MJ->at(1)) > TMath::Pi()){
	    dphijj_monoJet.back()->Fill(2*TMath::Pi()-fabs(jetphi_MJ->at(0)-jetphi_MJ->at(1)),weight);
	    dRjj_monoJet.back()->Fill(sqrt(TMath::Power(2*TMath::Pi()-fabs(jetphi_MJ->at(0)-jetphi_MJ->at(1)),2)+TMath::Power(jeteta_MJ->at(0)-jeteta_MJ->at(1),2)),weight);
	  }
	  else{
	    dphijj_monoJet.back()->Fill(fabs(jetphi_MJ->at(0)-jetphi_MJ->at(1)),weight);
	    dRjj_monoJet.back()->Fill(sqrt(TMath::Power(fabs(jetphi_MJ->at(0)-jetphi_MJ->at(1)),2)+TMath::Power(jeteta_MJ->at(0)-jeteta_MJ->at(1),2)),weight);
	  }	
	}
	
	if(jetphi_MJ->size() < 2){
	  mindphijj_monoJet.back()->Fill(mindphijj_monoJet.back()->GetBinCenter(1),weight);
	  mindphij1j_monoJet.back()->Fill(mindphij1j_monoJet.back()->GetBinCenter(1),weight);
	  mindRj1j_monoJet.back()->Fill(dRjj_monoJet.back()->GetBinCenter(1),weight);
	}
	else{
	  float minDphijj  = 2*TMath::Pi();
	  float minDphij1j = 2*TMath::Pi();
	  bool  isfound_jj  = false;
	  bool  isfound_j1j = false;
	  for(size_t ijet = 0 ; ijet < jetphi_MJ->size(); ijet++){
	    for(size_t jjet = ijet+1 ; jjet < jetphi_MJ->size(); jjet++){
	      float deltaPhi = fabs(jetphi_MJ->at(ijet)-jetphi_MJ->at(jjet));
	      if(jetpt_MJ->at(ijet) < 30 or jetpt_MJ->at(jjet) < 30) continue;
	      if(deltaPhi > TMath::Pi())
		deltaPhi = 2*TMath::Pi() - deltaPhi;
	      if(deltaPhi > 0 and deltaPhi < minDphijj){
		minDphijj  = deltaPhi;
		isfound_jj = true;
	      }
	    }
	  }
	  for(size_t jjet = 1 ; jjet < jetphi_MJ->size(); jjet++){
	    if(jetpt_MJ->at(0) < 30 or jetpt_MJ->at(jjet) < 30) continue;
	    float deltaPhi = fabs(jetphi_MJ->at(0)-jetphi_MJ->at(jjet));
	    if(deltaPhi > TMath::Pi())
	      deltaPhi = 2*TMath::Pi() - deltaPhi;
	    if(deltaPhi > 0 and deltaPhi < minDphij1j){
	      minDphij1j  = deltaPhi;
	      isfound_j1j = true;
	    }
	  }
	  
	  if(isfound_jj)
	    mindphijj_monoJet.back()->Fill(minDphijj,weight);
	  else 
	    mindphijj_monoJet.back()->Fill(mindphijj_monoJet.back()->GetBinCenter(1),weight);
	  
	  if(isfound_j1j){
	    mindphij1j_monoJet.back()->Fill(minDphij1j,weight);
	    mindRj1j_monoJet.back()->Fill(sqrt(minDphij1j*minDphij1j+TMath::Power(jeteta_MJ->at(0)-jeteta_MJ->at(1),2)),weight);
	  }
	  else{
	    mindphij1j_monoJet.back()->Fill(mindphij1j_monoJet.back()->GetBinCenter(1),weight);
	    mindRj1j_monoJet.back()->Fill(mindphij1j_monoJet.back()->GetBinCenter(1),weight);
	  }
	}
      }
      std::cout<<std::endl;
    }    
  }
    
  cout<<"Backgroud analysis"<<endl;
  TTreeReader myReader(backgroundZnunu);
  TTreeFormula* formula_bkg = new TTreeFormula("Znunu_formula",cut.c_str(),backgroundZnunu);  
  
  TH1F* boostedJetPt_bkg = NULL;
  TH1F* leadingJetEta_bkg = new TH1F("leadingJetEta_bkg","Z #rightarrow #nu #nu",20,-2.5,2.5);  
  TH1F* boostedJetEta_bkg = new TH1F("boostedJetEta_bkg","Z #rightarrow #nu #nu",20,-2.5,2.5);
  TH1F* boostedJetTau2Tau1_bkg = NULL;
  TH1F* boostedPrunedJetm_bkg = NULL;
  TH1F* leadingJetPt_bkg  = NULL;
  TH1F* leadingJetQGL_bkg = NULL;
  TH1F* njets_bkg = NULL;
  TH1F* njets_fwd_bkg = NULL;
  TH1F* ht_bkg    = NULL;
  TH1F* mT_bkg    = NULL;
  TH1F* met_bkg   =  NULL;
  TH1F* dphijj_bkg   = NULL;
  TH1F* mindphijj_bkg   = NULL;
  TH1F* mindphij1j_bkg   = NULL;
  TH1F* dRjj_bkg = NULL;
  TH1F* mindRj1j_bkg = NULL;

  if(category == 2){
    boostedJetPt_bkg = new TH1F("boostedJetPt_bkg","Z #rightarrow #nu #nu",bins_monoV_jetPt.size()-1,&bins_monoV_jetPt[0]);
    boostedJetTau2Tau1_bkg = new TH1F("boostedJetTau2Tau1_bkg","Z #rightarrow #nu #nu",bins_monoV_tau2tau1.size()-1,&bins_monoV_tau2tau1[0]);
    boostedPrunedJetm_bkg = new TH1F("boostedPrunedJetm_bkg","Z #rightarrow #nu #nu",bins_monoV_mpr.size()-1,&bins_monoV_mpr[0]);
    leadingJetPt_bkg  = new TH1F("leadingJetPt_bkg","Z #rightarrow #nu #nu",bins_monoV_jetPt.size()-1,&bins_monoV_jetPt[0]);  
    leadingJetQGL_bkg = new TH1F("leadingJetQGL_bkg","Z #rightarrow #nu #nu",bins_monoV_QGL.size()-1,&bins_monoV_QGL[0]);  
    njets_bkg = new TH1F("njets_bkg","Z #rightarrow #nu #nu",bins_monoV_njet.size()-1,&bins_monoV_njet[0]);
    njets_fwd_bkg = new TH1F("njets_fwd_bkg","Z #rightarrow #nu #nu",bins_monoV_njet.size()-1,&bins_monoV_njet[0]);
    ht_bkg = new TH1F("ht_bkg","Z #rightarrow #nu #nu",bins_monoV_HT.size()-1,&bins_monoV_HT[0]);
    mT_bkg = new TH1F("mT_bkg","Z #rightarrow #nu #nu",bins_monoV_mT.size()-1,&bins_monoV_mT[0]);
    met_bkg = new TH1F("met_bkg","Z #rightarrow #nu #nu",bins_monoV_met.size()-1,&bins_monoV_met[0]);
    dphijj_bkg = new TH1F("dphijj_bkg","Z #rightarrow #nu #nu",bins_monoV_dphiJJ.size()-1,&bins_monoV_dphiJJ[0]);
    mindphijj_bkg = new TH1F("mindphijj_bkg","Z #rightarrow #nu #nu",bins_monoV_dphiJJ.size()-1,&bins_monoV_dphiJJ[0]);
    mindphij1j_bkg = new TH1F("mindphij1j_bkg","Z #rightarrow #nu #nu",bins_monoV_dphiJJ.size()-1,&bins_monoV_dphiJJ[0]);
    dRjj_bkg = new TH1F("dRjj_bkg","Z #rightarrow #nu #nu",bins_monoV_dRJJ.size()-1,&bins_monoV_dRJJ[0]);
    mindRj1j_bkg = new TH1F("mindRj1j_bkg","Z #rightarrow #nu #nu",bins_monoV_dRJJ.size()-1,&bins_monoV_dRJJ[0]);

  }
  else{
    boostedJetPt_bkg = new TH1F("boostedJetPt_bkg","Z #rightarrow #nu #nu",bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]);
    boostedJetTau2Tau1_bkg = new TH1F("boostedJetTau2Tau1_bkg","Z #rightarrow #nu #nu",bins_monoJ_tau2tau1.size()-1,&bins_monoJ_tau2tau1[0]);
    boostedPrunedJetm_bkg = new TH1F("boostedPrunedJetm_bkg","Z #rightarrow #nu #nu",bins_monoJ_mpr.size()-1,&bins_monoJ_mpr[0]);
    leadingJetPt_bkg  = new TH1F("leadingJetPt_bkg","Z #rightarrow #nu #nu",bins_monoJ_jetPt.size()-1,&bins_monoJ_jetPt[0]);  
    leadingJetQGL_bkg = new TH1F("leadingJetQGL_bkg","Z #rightarrow #nu #nu",bins_monoJ_QGL.size()-1,&bins_monoJ_QGL[0]);  
    njets_bkg = new TH1F("njets_bkg","Z #rightarrow #nu #nu",bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]);
    njets_fwd_bkg = new TH1F("njets_fwd_bkg","Z #rightarrow #nu #nu",bins_monoJ_njet.size()-1,&bins_monoJ_njet[0]);
    ht_bkg = new TH1F("ht_bkg","Z #rightarrow #nu #nu",bins_monoJ_HT.size()-1,&bins_monoJ_HT[0]);
    mT_bkg = new TH1F("mT_bkg","Z #rightarrow #nu #nu",bins_monoJ_mT.size()-1,&bins_monoJ_mT[0]);
    met_bkg = new TH1F("met_bkg","Z #rightarrow #nu #nu",bins_monoJ_met.size()-1,&bins_monoJ_met[0]);
    dphijj_bkg = new TH1F("dphijj_bkg","Z #rightarrow #nu #nu",bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]);
    mindphijj_bkg = new TH1F("mindphijj_bkg","Z #rightarrow #nu #nu",bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]);
    mindphij1j_bkg = new TH1F("mindphij1j_bkg","Z #rightarrow #nu #nu",bins_monoJ_dphiJJ.size()-1,&bins_monoJ_dphiJJ[0]);    
    dRjj_bkg = new TH1F("dRjj_bkg","Z #rightarrow #nu #nu",bins_monoJ_dRJJ.size()-1,&bins_monoJ_dRJJ[0]);
    mindRj1j_bkg = new TH1F("mindRj1j_bkg","Z #rightarrow #nu #nu",bins_monoJ_dRJJ.size()-1,&bins_monoJ_dRJJ[0]);
  }

  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
  TTreeReaderValue<vector<double> > boostedJetm_pr  (myReader,"prunedJetm");
  
  TTreeReaderValue<vector<double> > jetpt    (myReader,"combinejetpt");
  TTreeReaderValue<vector<double> > jeteta   (myReader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetphi   (myReader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetQGL   (myReader,"combinejetQGL");
  TTreeReaderValue<unsigned int> njets  (myReader,"njets");
  TTreeReaderValue<unsigned int> nbjets  (myReader,"nbjetslowpt");
  TTreeReaderValue<vector<double> > fwdjetpt  (myReader,"forwardjetpt");
  TTreeReaderValue<double> ht  (myReader,"ht");
  TTreeReaderValue<double> met  (myReader,"t1pfmet");
  TTreeReaderValue<double> metphi  (myReader,"t1pfmetphi");    

  // weight for xsec to account for binning in HT of the sample
  TTreeReaderValue<double> xsec         (myReader,"xsec");
  TTreeReaderValue<double> wgt          (myReader,"wgt");
  TTreeReaderValue<double> wgtbtag      (myReader,"wgtbtag");
  TTreeReaderValue<double> wgtsum       (myReader,"wgtsum");
  TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");
  TTreeReaderValue<double> wzpt  (myReader,"wzpt");

  // get k-factors NLO                                                                                                                                                           
  TFile kffile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_Wseparated_24bins.root");
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

  vector<TH1*> khists;
  khists.push_back(znlohist);
  khists.push_back(zewkhist);
  
  int iEvent = 0;
  cout<<"Run on Znunu events "<<backgroundZnunu->GetEntries()<<endl;
  while(myReader.Next()){
    backgroundZnunu->GetEntry(iEvent);
    cout.flush();
    if(iEvent %25000 == 0) cout<<"\r"<<"iEvent "<<100*float(iEvent)/backgroundZnunu->GetEntries()<<" %";
    iEvent++;
    if(not formula_bkg->EvalInstance()) continue;

    Double_t kwgt = 1.0;
    for (unsigned i = 0; i < khists.size(); i++) {
      if(*wzpt < khists[i]->GetBinLowEdge(1))
	*wzpt  = khists[i]->GetBinLowEdge(1)+1.;
      else if(*wzpt > khists[i]->GetBinLowEdge(khists[i]->GetNbinsX()+1))
	*wzpt = khists[i]->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1.;

      kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(*wzpt));
    }

    float weight = (*xsec)*(lumi)*(*wgt)*(*wgtbtag)*(puhist->GetBinContent(*nvtx))*(trmhist->GetBinContent(trmhist->FindBin(*met)))*kwgt/(*wgtsum);

    if(boostedJetpt->size()>0)
      boostedJetPt_bkg->Fill(boostedJetpt->at(0),weight);
    if(boostedJeteta->size()>0)
      boostedJetEta_bkg->Fill(boostedJeteta->at(0),weight);
    if(boostedJettau2->size()>0 and boostedJettau1->size()>0 )
      boostedJetTau2Tau1_bkg->Fill(boostedJettau2->at(0)/boostedJettau1->at(0),weight);
    if(boostedJetm_pr->size()>0)
      boostedPrunedJetm_bkg->Fill(boostedJetm_pr->at(0),weight);
    if(jetpt->size()>0)
      leadingJetPt_bkg->Fill(jetpt->at(0),weight);
    if(jeteta->size()>0)
      leadingJetEta_bkg->Fill(jeteta->at(0),weight);
    if(jetQGL->size()>0)
      leadingJetQGL_bkg->Fill(jetQGL->at(0),weight);
    
    njets_bkg->Fill(*njets,weight);
    njets_fwd_bkg->Fill(*njets+fwdjetpt->size(),weight);
    ht_bkg->Fill(*ht,weight);
    met_bkg->Fill(*met,weight);
    
    if(jetpt->size()>0){
      float deltaPhi = fabs(jetphi->at(0)-*metphi);
      if(deltaPhi > TMath::Pi())
	deltaPhi = fabs(2*TMath::Pi() - deltaPhi);
      mT_bkg->Fill(sqrt(2*jetpt->at(0)*(*met)*(1-cos(deltaPhi))),weight);
    }
   

    //////
    if(jetphi->size() < 2){
      dphijj_bkg->Fill(dphijj_bkg->GetBinCenter(1),weight);
      dRjj_bkg->Fill(dRjj_bkg->GetBinCenter(1),weight);
    }
    else{
      if(fabs(jetphi->at(0)-jetphi->at(1)) > TMath::Pi()){
	dphijj_bkg->Fill(2*TMath::Pi()-fabs(jetphi->at(0)-jetphi->at(1)),weight);
	dRjj_bkg->Fill(sqrt(TMath::Power(2*TMath::Pi()-fabs(jetphi->at(0)-jetphi->at(1)),2)+TMath::Power(jeteta->at(0)-jeteta->at(1),2)),weight);
      }
      else{
	dphijj_bkg->Fill(fabs(jetphi->at(0)-jetphi->at(1)),weight);
	dRjj_bkg->Fill(sqrt(TMath::Power(fabs(jetphi->at(0)-jetphi->at(1)),2)+TMath::Power(jeteta->at(0)-jeteta->at(1),2)),weight);
      }	
    }
    
    if(jetphi->size() < 2){
      mindphijj_bkg->Fill(mindphijj_bkg->GetBinCenter(1),weight);
      mindphij1j_bkg->Fill(mindphij1j_bkg->GetBinCenter(1),weight);
      mindRj1j_bkg->Fill(dRjj_bkg->GetBinCenter(1),weight);
    }
    else{
      float minDphijj  = 2*TMath::Pi();
      float minDphij1j = 2*TMath::Pi();
      bool  isfound_jj  = false;
      bool  isfound_j1j = false;
      for(size_t ijet = 0 ; ijet < jetphi->size(); ijet++){
	for(size_t jjet = ijet+1 ; jjet < jetphi->size(); jjet++){
	  if(jetpt->at(ijet) < 30 or jetpt->at(jjet) < 30) continue;
	  float deltaPhi = fabs(jetphi->at(ijet)-jetphi->at(jjet));
	    if(deltaPhi > TMath::Pi())
	      deltaPhi = 2*TMath::Pi() - deltaPhi;
	    if(deltaPhi > 0 and deltaPhi < minDphijj){
	      minDphijj  = deltaPhi;
	      isfound_jj = true;
	    }
	}
      }
      for(size_t jjet = 1 ; jjet < jetphi->size(); jjet++){
	if(jetpt->at(0) < 30 or jetpt->at(jjet) < 30) continue;
	float deltaPhi = fabs(jetphi->at(0)-jetphi->at(jjet));
	if(deltaPhi > TMath::Pi())
	  deltaPhi = 2*TMath::Pi() - deltaPhi;
	if(deltaPhi > 0 and deltaPhi < minDphij1j){
	  minDphij1j  = deltaPhi;
	  isfound_j1j = true;
	}
      }
      
      if(isfound_jj)
	mindphijj_bkg->Fill(minDphijj,weight);
      else 
	mindphijj_bkg->Fill(mindphijj_bkg->GetBinCenter(1),weight);
      
      if(isfound_j1j){
	mindphij1j_bkg->Fill(minDphij1j,weight);
	mindRj1j_bkg->Fill(sqrt(minDphij1j*minDphij1j+TMath::Power(jeteta->at(0)-jeteta->at(1),2)),weight);
      }
      else{
	  mindphij1j_bkg->Fill(mindphij1j_bkg->GetBinCenter(1),weight);
	  mindRj1j_bkg->Fill(mindphij1j_bkg->GetBinCenter(1),weight);
      }
    }
  }
  std::cout<<std::endl;
  // make single plots
  if(category == 0 and not displayOnlyMonoV){
    
    makeShapePlots(cCanvas,boostedJetPt_monoV,boostedJetPt_monoJet,boostedJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T}^{AK8} (GeV)");
    makeShapePlots(cCanvas,boostedJetEta_monoV,boostedJetEta_monoJet,boostedJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta^{AK8}");    
    makeShapePlots(cCanvas,boostedJetTau2Tau1_monoV,boostedJetTau2Tau1_monoJet,boostedJetTau2Tau1_bkg,
		   legend,outputDirectory,"jetTau2Tau1","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,boostedPrunedJetm_monoV,boostedPrunedJetm_monoJet,boostedPrunedJetm_bkg,
		   legend,outputDirectory,"mpruned","m_{pruned} (GeV)");
    
    makeShapePlots(cCanvas,leadingJetPt_monoV,leadingJetPt_monoJet,leadingJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T} (GeV)");
    makeShapePlots(cCanvas,leadingJetEta_monoV,leadingJetEta_monoJet,leadingJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta");    
    makeShapePlots(cCanvas,leadingJetQGL_monoV,leadingJetQGL_monoJet,leadingJetQGL_bkg,
		   legend,outputDirectory,"jetQGL","QGL");    
    makeShapePlots(cCanvas,njets_monoV,njets_monoJet,njets_bkg,
		   legend,outputDirectory,"njet","N_{jet}");    
    makeShapePlots(cCanvas,njets_fwd_monoV,njets_fwd_monoJet,njets_fwd_bkg,
		   legend,outputDirectory,"njet_fwd","N_{jet}");    
    makeShapePlots(cCanvas,ht_monoV,ht_monoJet,ht_bkg,
		 legend,outputDirectory,"ht","H_{T} (GeV)");    
    makeShapePlots(cCanvas,mT_monoV,mT_monoJet,mT_bkg,
		   legend,outputDirectory,"mT","m_{T} (GeV)");    
    makeShapePlots(cCanvas,met_monoV,met_monoJet,met_bkg,
		   legend,outputDirectory,"met","E_{T}^{miss} (GeV)");    
    makeShapePlots(cCanvas,dphijj_monoV,dphijj_monoJet,dphijj_bkg,
		   legend,outputDirectory,"dphijj","#Delta#phi(j_{1},j_{2})");    
    makeShapePlots(cCanvas,mindphijj_monoV,mindphijj_monoJet,mindphijj_bkg,
		   legend,outputDirectory,"mindphijj","min(#Delta#phi(j,j))");    
    makeShapePlots(cCanvas,mindphij1j_monoV,mindphij1j_monoJet,mindphij1j_bkg,
		   legend,outputDirectory,"mindphij1j","min(#Delta#phi(j_{1},j))");    

    makeShapePlots(cCanvas,dRjj_monoV,dRjj_monoJet,dRjj_bkg,
		   legend,outputDirectory,"dRjj","#DeltaR(j_{1},j_{2})");    
    makeShapePlots(cCanvas,mindRj1j_monoV,mindRj1j_monoJet,mindRj1j_bkg,
		   legend,outputDirectory,"mindRj1j","min(#DeltaR(j_{1},j))");    
  }
  else if(category == 1){

    makeShapePlots(cCanvas,boostedJetPt_monoJet,boostedJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T}^{AK8} (GeV)");
    makeShapePlots(cCanvas,boostedJetEta_monoJet,boostedJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta^{AK8}");    
    makeShapePlots(cCanvas,boostedJetTau2Tau1_monoJet,boostedJetTau2Tau1_bkg,
		   legend,outputDirectory,"jetTau2Tau1","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,boostedPrunedJetm_monoJet,boostedPrunedJetm_bkg,
		   legend,outputDirectory,"mpruned","m_{pruned} (GeV)");
    
    makeShapePlots(cCanvas,leadingJetPt_monoJet,leadingJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T} (GeV)");
    makeShapePlots(cCanvas,leadingJetEta_monoJet,leadingJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta");    
    makeShapePlots(cCanvas,leadingJetQGL_monoJet,leadingJetQGL_bkg,
		   legend,outputDirectory,"jetQGL","QGL");    
    makeShapePlots(cCanvas,njets_monoJet,njets_bkg,
		   legend,outputDirectory,"njet","N_{jet}");    
    makeShapePlots(cCanvas,njets_fwd_monoJet,njets_fwd_bkg,
		   legend,outputDirectory,"njet_fwd","N_{jet}");    
    makeShapePlots(cCanvas,ht_monoJet,ht_bkg,
		 legend,outputDirectory,"ht","H_{T} (GeV)");    
    makeShapePlots(cCanvas,mT_monoJet,mT_bkg,
		   legend,outputDirectory,"mT","m_{T} (GeV)");    
    makeShapePlots(cCanvas,met_monoJet,met_bkg,
		   legend,outputDirectory,"met","E_{T}^{miss} (GeV)");    
    makeShapePlots(cCanvas,dphijj_monoJet,dphijj_bkg,
		   legend,outputDirectory,"dphijj","#Delta#phi(j_{1},j_{2})");    
    makeShapePlots(cCanvas,mindphijj_monoJet,mindphijj_bkg,
		   legend,outputDirectory,"mindphijj","min(#Delta#phi(j,j))");    
    makeShapePlots(cCanvas,mindphij1j_monoJet,mindphij1j_bkg,
		   legend,outputDirectory,"mindphij1j","min(#Delta#phi(j_{1},j))");    

    makeShapePlots(cCanvas,dRjj_monoJet,dRjj_bkg,
		   legend,outputDirectory,"dRjj","#DeltaR(j_{1},j_{2})");    
    makeShapePlots(cCanvas,mindRj1j_monoJet,mindRj1j_bkg,
		   legend,outputDirectory,"mindRj1j","min(#DeltaR(j_{1},j))");    

  }
  else if(category == 2 or (category == 0 and displayOnlyMonoV)){

    makeShapePlots(cCanvas,boostedJetPt_monoV,boostedJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T}^{AK8} (GeV)");
    makeShapePlots(cCanvas,boostedJetEta_monoV,boostedJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta^{AK8}");    
    makeShapePlots(cCanvas,boostedJetTau2Tau1_monoV,boostedJetTau2Tau1_bkg,
		   legend,outputDirectory,"jetTau2Tau1","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,boostedPrunedJetm_monoV,boostedPrunedJetm_bkg,
		   legend,outputDirectory,"mpruned","m_{pruned} (GeV)");
    
    makeShapePlots(cCanvas,leadingJetPt_monoV,leadingJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T} (GeV)");
    makeShapePlots(cCanvas,leadingJetEta_monoV,leadingJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta");    
    makeShapePlots(cCanvas,leadingJetQGL_monoV,leadingJetQGL_bkg,
		   legend,outputDirectory,"jetQGL","QGL");    
    makeShapePlots(cCanvas,njets_monoV,njets_bkg,
		   legend,outputDirectory,"njet","N_{jet}");    
    makeShapePlots(cCanvas,njets_fwd_monoV,njets_fwd_bkg,
		   legend,outputDirectory,"njet_fwd","N_{jet}");    
    makeShapePlots(cCanvas,ht_monoV,ht_bkg,
		 legend,outputDirectory,"ht","H_{T} (GeV)");    
    makeShapePlots(cCanvas,mT_monoV,mT_bkg,
		   legend,outputDirectory,"mT","m_{T} (GeV)");    
    makeShapePlots(cCanvas,met_monoV,met_bkg,
		   legend,outputDirectory,"met","E_{T}^{miss} (GeV)");    
    makeShapePlots(cCanvas,dphijj_monoV,dphijj_bkg,
		   legend,outputDirectory,"dphijj","#Delta#phi(j_{1},j_{2})");    
    makeShapePlots(cCanvas,mindphijj_monoV,mindphijj_bkg,
		   legend,outputDirectory,"mindphijj","min(#Delta#phi(j_{1},j_{2}))");    
    makeShapePlots(cCanvas,mindphij1j_monoV,mindphij1j_bkg,
		   legend,outputDirectory,"mindphij1j","min(#Delta#phi(j_{1},j))");    
    makeShapePlots(cCanvas,dRjj_monoV,dRjj_bkg,
		   legend,outputDirectory,"dRjj","#DeltaR(j_{1},j_{2})");    
    makeShapePlots(cCanvas,mindRj1j_monoV,mindRj1j_bkg,
		   legend,outputDirectory,"mindRj1j","min(#DeltaR(j_{1},j))");    

  }

  return;
}



void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList_MV, vector<TH1F*> histoList_MJ,
		    TH1F* histoBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTitle){

  legend->Clear();

  TLatex * tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.276,0.92,"Preliminary Simulation");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  float maxYScale = 0;
  float minYScale = 0;

  float integral = histoBkg->Integral();
  histoBkg->Scale(1./integral);
  maxYScale = histoBkg->GetMaximum();
  minYScale = histoBkg->GetMinimum();

  for(size_t iHisto = 0; iHisto < histoList_MV.size(); iHisto++){

    float integral = histoList_MV.at(iHisto)->Integral();
    histoList_MV.at(iHisto)->Scale(1./integral);

    if(histoList_MV.at(iHisto)->GetMaximum() > maxYScale)
      maxYScale = histoList_MV.at(iHisto)->GetMaximum();
    if(histoList_MV.at(iHisto)->GetMinimum() < minYScale)
      minYScale = histoList_MV.at(iHisto)->GetMinimum();

  }

  for(size_t iHisto = 0; iHisto < histoList_MJ.size(); iHisto++){

    float integral = histoList_MJ.at(iHisto)->Integral();
    histoList_MJ.at(iHisto)->Scale(1./integral);

    if(histoList_MJ.at(iHisto)->GetMaximum() > maxYScale)
      maxYScale = histoList_MJ.at(iHisto)->GetMaximum();

    if(histoList_MJ.at(iHisto)->GetMaximum() > minYScale)
      minYScale = histoList_MJ.at(iHisto)->GetMinimum();
  }

  histoBkg->SetLineColor(kBlack);
  histoBkg->SetLineWidth(2);
  histoBkg->SetFillColor(kGray);
  histoBkg->GetXaxis()->SetTitleOffset(1.15);
  histoBkg->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoBkg->GetYaxis()->SetRangeUser(0,maxYScale*1.5);
  histoBkg->GetYaxis()->SetTitle("a.u.");
  histoBkg->GetYaxis()->SetTitleOffset(1.15);
  histoBkg->Draw("hist");

  legend->AddEntry(histoBkg,histoBkg->GetTitle(),"FL");

  int icolor = 0;

  for(size_t iHisto = 0; iHisto < histoList_MV.size(); iHisto++){
    
    histoList_MV.at(iHisto)->SetLineWidth(2);
    histoList_MV.at(iHisto)->SetLineColor(color[icolor]);

    histoList_MV.at(iHisto)->SetMarkerStyle(20);
    histoList_MV.at(iHisto)->SetMarkerSize(1.5);
        
    
    legend->AddEntry(histoList_MV.at(iHisto),histoList_MV.at(iHisto)->GetTitle(),"l");
    
    histoList_MV.at(iHisto)->Draw("histsame");
    
  }

  for(size_t iHisto = 0; iHisto < histoList_MJ.size(); iHisto++){
    
    histoList_MJ.at(iHisto)->SetLineWidth(2);
    histoList_MJ.at(iHisto)->SetLineColor(color[icolor]);
    icolor++;

    histoList_MJ.at(iHisto)->SetMarkerStyle(20);
    histoList_MJ.at(iHisto)->SetMarkerSize(1.5);

    legend->AddEntry(histoList_MJ.at(iHisto),histoList_MJ.at(iHisto)->GetTitle(),"l");

    histoList_MJ.at(iHisto)->Draw("histsame");
    
  }

  legend->Draw("same");
  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  
  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB.png").c_str(),"png");
  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB.pdf").c_str(),"pdf");

  if(minYScale != 0)
    histoBkg->GetYaxis()->SetRangeUser(minYScale*0.5,maxYScale*500);
  else
    histoBkg->GetYaxis()->SetRangeUser(0.0001,maxYScale*500);
  cCanvas->SetLogy();

  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB_log.png").c_str(),"png");
  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB_log.pdf").c_str(),"pdf");
  cCanvas->SetLogy(0);  
}



void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList_MV,
		    TH1F* histoBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTitle){

  legend->Clear();

  TLatex * tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.286,0.92,"Preliminary Simulation");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  float maxYScale = 0;
  float minYScale = 0;

  float integral = histoBkg->Integral();
  histoBkg->Scale(1./integral);
  maxYScale = histoBkg->GetMaximum();
  minYScale = histoBkg->GetMinimum();

  for(size_t iHisto = 0; iHisto < histoList_MV.size(); iHisto++){

    float integral = histoList_MV.at(iHisto)->Integral();
    histoList_MV.at(iHisto)->Scale(1./integral);

    if(histoList_MV.at(iHisto)->GetMaximum() > maxYScale)
      maxYScale = histoList_MV.at(iHisto)->GetMaximum();
    if(histoList_MV.at(iHisto)->GetMinimum() < minYScale)
      minYScale = histoList_MV.at(iHisto)->GetMinimum();

  }


  histoBkg->SetLineColor(kBlack);
  histoBkg->SetLineWidth(2);
  histoBkg->SetFillColor(kGray);
  histoBkg->GetXaxis()->SetTitleOffset(0.97);
  histoBkg->GetXaxis()->SetNdivisions(509);
  histoBkg->GetYaxis()->SetNdivisions(509);
  histoBkg->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoBkg->GetYaxis()->SetRangeUser(0,maxYScale*1.5);
  histoBkg->GetYaxis()->SetTitle("a.u.");
  histoBkg->GetXaxis()->SetTitleSize(0.05);
  histoBkg->GetYaxis()->SetTitleSize(0.05);
  histoBkg->GetYaxis()->SetTitleOffset(1.27);
  histoBkg->Draw("hist");

  legend->AddEntry(histoBkg,histoBkg->GetTitle(),"F");

  int icolor = 0;

  for(size_t iHisto = 0; iHisto < histoList_MV.size(); iHisto++){
    
    histoList_MV.at(iHisto)->SetLineWidth(2);
    histoList_MV.at(iHisto)->SetLineColor(color[icolor]);
    icolor++;

    histoList_MV.at(iHisto)->SetMarkerStyle(20);
    histoList_MV.at(iHisto)->SetMarkerSize(1.5);

    string legendName;
    TString title (histoList_MV.at(iHisto)->GetTitle());
    if(title.Contains("Vector"))
      legendName += "Mono-W Vector, ";
    else if(title.Contains("Axial"))
      legendName += "Mono-W Axial, ";
    else if(title.Contains("Scalar"))
      legendName += "Mono-W Scalar, ";
    else if(title.Contains("Pseudoscalar"))
      legendName += "Mono-W Pseudoscalar, ";

    std::stringstream test(title.Data());
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(test, segment, '_')){
      seglist.push_back(segment);
    }

    for(auto seg : seglist){
      if(TString(seg).Contains("Mphi-"))
	legendName += "m_{med} = "+TString(seg).ReplaceAll("Mphi-","")+" GeV, ";
      else if(TString(seg).Contains("Mchi-"))
	legendName += "m_{DM} = "+TString(seg).ReplaceAll("Mchi-","")+" GeV";
    }
    
    legend->AddEntry(histoList_MV.at(iHisto),legendName.c_str(),"l");
    histoList_MV.at(iHisto)->Draw("histsame");
    
  }


  legend->Draw("same");
  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  
  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB.pdf").c_str(),"pdf");
  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB.png").c_str(),"png");

  if(minYScale != 0)
    histoBkg->GetYaxis()->SetRangeUser(minYScale*0.5,maxYScale*1.5);
  else
    histoBkg->GetYaxis()->SetRangeUser(0.0001,maxYScale*500);
  cCanvas->SetLogy();

  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs((outputDirectory+"/"+plotName+"_SandB_log.png").c_str(),"png");
  cCanvas->SetLogy(0);  
}
