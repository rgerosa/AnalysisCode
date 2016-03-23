#ifndef HISTOUTILS_H
#define HISTOUTILS_H

#include <vector>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TString.h"

using namespace std;

class signalSample{
  
 public:
  signalSample(string a,string b, string c){
    interaction = a;
    mediatorMass = b;
    dmMass = c;
  }
  

  string interaction;
  string mediatorMass;
  string dmMass;
};

class VectorSorter{
 public:
  bool operator ()(const TLorentzVector & i, const TLorentzVector & j) const {
    return (i.Pt() > j.Pt());
  }
};

// define binnings for the different observables                                 
vector<double> bins_monoV_met         = {250.,300.,350.,400.,500.,600.,750.,1000.};
vector<double> bins_monoV_met_v2      = {250.,300.,350.,400.,500.,600.,700.,1000.};
vector<double> bins_monoV_met_v3      = {250.,300.,350.,400.,500.,600.,800.,1000.};
vector<double> bins_substructure_met  = {250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_monoJ_met         = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
vector<double> bins_monoJ_met_v2      = {200.,250.,300.,350.,400.,500.,600.,1000.};

vector<double> bins_monoJ_dphiJJ      = {0.,0.25,0.5,0.65,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.14};
vector<double> bins_monoV_dphiJJ      = {0.,0.25,0.5,0.65,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.14};

vector<double> bins_monoV_mT          = {50.,100.,150.,200.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_substructure_mT   = {50.,100.,150.,200.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_monoJ_mT          = {50.,80.,110.,140.,170,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1350,1550,1750,2000};

vector<double> bins_monoV_mpr         = {65.,67.5,70.,72.5,75.,77.5,80.,82.5,85.,87.5,90.,92.5,95.,97.5,100.,102.5,105.};
vector<double> bins_monoV_mpr_v2      = {65.,73.,81.,89.,97.,105.};

vector<double> bins_monoJ_mpr         = {0.,3.,6.,9.,12.,15.,18.,21.,24.,27.,30.,33.,36.,39.,42.,45.,48.,51.,54.,57.,60.,64.,68.,72.,76.,80.,84.,88.,92.,96.,100.};
vector<double> bins_substructure_mpr  = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};

vector<double> bins_monoV_njet        = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
vector<double> bins_monoJ_njet        = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
vector<double> bins_substructure_njet = {0.,1.,2.,3.,4.,5.,6.,7.,8.};

vector<double> bins_monoV_HT          = {0,50.,100.,150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
vector<double> bins_monoJ_HT          = {0,50.,100.,150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
vector<double> bins_substructure_HT   = {0,50.,100.,150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};

vector<double> bins_monoV_jetPt        = {200.,225.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_monoJ_jetPt        = {100.,120.,140.,160.,180.,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
vector<double> bins_substructure_jetPt = {100.,120.,140.,160.,180.,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};

vector<double> bins_monoV_bosonPt    = {200.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_monoV_bosonPt_v2 = {50.,70.,90.,120.,150.,180.,210.,230.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_monoJ_bosonPt    = {50.,70.,90.,120.,150.,180.,210.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250.};
vector<double> bins_monoJ_bosonPt_v2 = {200.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_substructure_bosonPt = {50.,70.,90.,120.,150.,180.,210.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250.};

vector<double> bins_monoV_QGL        = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};
vector<double> bins_monoJ_QGL        = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};
vector<double> bins_substructure_QGL = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};

vector<double> bins_monoV_tau2tau1        = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.90,1.};
vector<double> bins_monoJ_tau2tau1        = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.90,1.};
vector<double> bins_substructure_tau2tau1 = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.90,1.};

vector<double> bins_monoV_nvtx        = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};
vector<double> bins_monoJ_nvtx        = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};
vector<double> bins_substructure_nvtx = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};

vector<double> bins_monoV_chfrac        = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
vector<double> bins_monoJ_chfrac        = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
vector<double> bins_substructure_chfrac = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};

vector<double> bins_monoV_nhfrac        = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
vector<double> bins_monoJ_nhfrac        = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
vector<double> bins_substructure_nhfrac = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};

vector<double> bins_monoJ_btagCSV = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};
vector<double> bins_monoV_btagCSV = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};

vector<double> bins_monoJ_met_2D      = {200.,250.,300.,350,400.,500.,600.,750.,900.,1050.,1300.};
vector<double> bins_monoJ_mpruned_2D  = {0.,5.,10.,15.,20.,25.,30.,35.,45.,55.,65.,75.,85.,95.,105.};
vector<double> bins_monoJ_tau2tau1_2D = {0.,0.15,0.3,0.4,0.5,0.7,0.8,0.9,1.};
vector<double> bins_monoJ_njet_2D     = {1,2,3,10};
vector<double> bins_monoJ_ht_2D       = {50.,300.,650.,950.,2000.};
vector<double> bins_monoJ_mT_2D       = {50.,450.,650.,950.,2000.};
vector<double> bins_monoJ_dphiJJ_2D   = {0.,0.25,1.3,3.14};
vector<double> bins_monoJ_QGL_2D      = {0.,0.15,0.50,0.85,1.};

// binning selections                                                                                                                                                          
vector<double> selectBinning (string observable, int category){

  if(observable == "met" and category <= 1 and category >= 0)
    return bins_monoJ_met;
  else if(observable == "met_v2" and category <= 1 and category >= 0)
    return bins_monoJ_met_v2;
  else if(observable == "met" and category > 1 and category <=3)
    return bins_monoV_met;
  else if(observable == "met_v2" and category > 1 and category <=3)
    return bins_monoV_met_v2;
  else if(observable == "met_v3" and category > 1 and category <=3)
    return bins_monoV_met_v3;
  else if(observable == "met" and category > 3)
    return bins_substructure_met;

  if(observable == "mT" and category <= 1)
    return bins_monoJ_mT;
  else if(observable == "mT" and category > 1 and category <=3)
    return bins_monoV_mT;
  else if(observable == "mT" and category > 3)
    return bins_substructure_mT;

  else if((observable == "mpruned" or observable == "msoftdrop" or observable == "mraw") and category <=1)
    return bins_monoJ_mpr;
  else if((observable == "mpruned" or observable == "msoftdrop" or observable == "mraw") and category > 1 and category <=3)
    return bins_monoV_mpr;
  else if((observable == "mpruned" or observable == "msoftdrop" or observable == "mraw") and category > 3)
    return bins_substructure_mpr;
  else if((observable == "mpruned_v2" or observable == "msoftdrop_v2" or observable == "mraw_v2") and category == 2)
    return bins_monoV_mpr_v2;

  else if(observable == "tau2tau1" and category <= 1)
    return bins_monoJ_tau2tau1;
  else if(observable == "tau2tau1" and category > 1 and category <=3)
    return bins_monoV_tau2tau1;
  else if(observable == "tau2tau1" and category > 1)
    return bins_substructure_tau2tau1;

  else if(observable == "njet" and category <= 1)
    return bins_monoJ_njet;
  else if(observable == "njet" and category > 1 and category <= 3)
    return bins_monoV_njet;
  else if(observable == "njet" and category > 3)
    return bins_substructure_njet;

  else if((observable == "nbjet" or observable == "nbjet_hpt" or observable == "nbjet_hpt_loose") and category <= 1)
    return bins_monoJ_njet;
  else if((observable == "nbjet" or observable == "nbjet_hpt" or observable == "nbjet_hpt_loose") and category > 1 and category <= 3)
    return bins_monoV_njet;
  else if((observable == "nbjet" or observable == "nbjet_hpt" or observable == "nbjet_hpt_loose") and category > 3)
    return bins_substructure_njet;
  
  else if(observable == "ht" and category <= 1)
    return bins_monoJ_HT;
  else if(observable == "ht" and category > 1 and category <= 3)
    return bins_monoV_HT;
  else if(observable == "ht" and category > 3)
    return bins_substructure_HT;

  else if(observable == "nvtx" and category <= 1)
    return bins_monoJ_nvtx;
  else if(observable == "nvtx" and category > 1 and category <= 3)
    return bins_monoV_nvtx;
  else if(observable == "nvtx" and category > 3)
    return bins_substructure_nvtx;

  else if(observable == "chfrac" and category <= 1)
    return bins_monoJ_chfrac;
  else if(observable == "chfrac" and category > 1 and category <= 3)
    return bins_monoV_chfrac;
  else if(observable == "chfrac" and category > 3)
    return bins_substructure_chfrac;

  else if(observable == "nhfrac" and category <= 1)
    return bins_monoJ_nhfrac;
  else if(observable == "nhfrac" and category > 1 and category <= 3)
    return bins_monoV_nhfrac;
  else if(observable == "nhfrac" and category > 3)
    return bins_substructure_nhfrac;

  else if(observable == "bosonPt" and category <= 1)
    return bins_monoJ_bosonPt;
  else if(observable == "bosonPt_v2" and category <= 1)
    return bins_monoJ_bosonPt_v2;
  else if(observable == "bosonPt" and category > 1 and category <= 3)
    return bins_monoV_bosonPt;
  else if(observable == "bosonPt" and category > 3)
    return bins_substructure_bosonPt;

  else if((observable == "jetPt" or observable == "boostedJetPt") and category <= 1)
    return bins_monoJ_jetPt;
  else if((observable == "jetPt" or observable == "boostedJetPt") and category > 1 and category <= 3)
    return bins_monoV_jetPt;
  else if((observable == "jetPt" or observable == "boostedJetPt") and category > 3)
    return bins_substructure_jetPt;

  else if((TString(observable).Contains("btag") or TString(observable).Contains("CSV")) and category <=1)
    return bins_monoJ_btagCSV;
  else if((TString(observable).Contains("btag") or TString(observable).Contains("CSV")) and category > 1)
    return bins_monoV_btagCSV;
  
  else if(TString(observable).Contains("QGL") and category <= 1)
    return bins_monoJ_QGL;
  else if(TString(observable).Contains("QGL") and category > 1 and category <= 3)
    return bins_monoV_QGL;
  else if(TString(observable).Contains("QGL") and category > 3)
    return bins_monoV_QGL;

  else if(observable == "dphiJJ" and category <= 1)
    return bins_monoJ_dphiJJ;
  else if(observable == "dphiJJ" and category > 1)
    return bins_monoV_dphiJJ;


  vector<double> dummy;
  return dummy;

}

struct bin2D {
  
  vector<double> binX;
  vector<double> binY;

};

bin2D selectBinning2D (string observable, int category){

  bin2D bins;

  if(observable == "met_mpruned" and category <= 1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_mpruned_2D;
    return bins;
  }
  else if(observable == "met_tau2tau1" and category <= 1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_tau2tau1_2D;
    return bins;
  }
  else if(observable == "mpruned_tau2tau1" and category <=1){
    bins.binX = bins_monoJ_mpruned_2D;
    bins.binY = bins_monoJ_tau2tau1_2D;
    return bins;
  }
  else if(observable == "met_njet" and category <=1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_njet_2D;
    return bins;
  }
  else if(observable == "met_ht" and category <=1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_ht_2D;
    return bins;
  }
  else if(observable == "met_mT" and category <=1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_mT_2D;
    return bins;
  }

  else if(observable == "met_QGL" and category <=1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_QGL_2D;
    return bins;    
  }

  else if(observable == "met_dphiJJ" and category <=1){
    bins.binX = bins_monoJ_met_2D;
    bins.binY = bins_monoJ_dphiJJ_2D;
    return bins;    
  }
  
  return bins;
}

// to smooth empty bins in TH1
void smoothEmptyBins(TH1* hist, int nsteps = 2){

  for(int iBin = 1 ; iBin <= hist->GetNbinsX(); iBin++){
    if(hist->GetBinContent(iBin) <= 0 or hist->GetBinError(iBin) >= hist->GetBinContent(iBin)){
      float average = 0.;
      for(int jBin = iBin -nsteps; jBin < iBin+nsteps; jBin++){
        if(jBin == iBin) continue;
        if(jBin > 0 and jBin <= hist->GetNbinsX()){
          average += hist->GetBinContent(jBin);
        }
      }
      hist->SetBinContent(iBin,average/(nsteps*2));
      hist->SetBinError(iBin,hist->GetBinContent(iBin)/(nsteps*2));
    }
  }
}

// to smooth empty bins in TH2
void smoothEmptyBins(TH2* hist, int nsteps = 1){
  for(int xBin = 1 ; xBin <= hist->GetNbinsX(); xBin++){
    for(int yBin = 1 ; yBin <= hist->GetNbinsY(); yBin++){
      if(hist->GetBinContent(xBin,yBin) <= 0 or hist->GetBinError(xBin,yBin) >= hist->GetBinContent(xBin,yBin)){
        float average = 0.;
        for(int jBin = xBin -nsteps; jBin < xBin+nsteps; jBin++){
          for(int kBin = yBin -nsteps; kBin < yBin+nsteps; kBin++){
            if(jBin == xBin and kBin == yBin) continue;
            if((jBin > 0 and jBin <= hist->GetNbinsX()) and (kBin > 0 and kBin <= hist->GetNbinsY())){
              average += hist->GetBinContent(jBin,kBin);
            }
          }
        }
        hist->SetBinContent(xBin,yBin,average/((nsteps*2+1)*(nsteps*2+1)-1));
        hist->SetBinError(xBin,yBin,hist->GetBinContent(xBin,yBin)*2);
      }
    }
  }
}

// make the average of two histograms (TH1)
void makeAverage(TH1* histo, TH1* histo_2){

  if(histo_2->Integral() != 0){ // sanity check to understand if the second has been filled
    for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
      histo->SetBinContent(iBin+1,(histo->GetBinContent(iBin+1)+histo_2->GetBinContent(iBin+1))/2);
      histo->SetBinError(iBin+1,sqrt(histo->GetBinError(iBin+1)*histo->GetBinError(iBin+1)+histo_2->GetBinError(iBin+1)*histo_2->GetBinError(iBin+1))*0.5);
    }
  }
}


// make the average of two histograms (TH2)
void makeAverage(TH2* histo, TH2* histo_2){

  if(histo_2->Integral() != 0){ // sanity check to understand if the second has been filled
    for(int iBinX = 0; iBinX < histo->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histo->GetNbinsY(); iBinY++){
	histo->SetBinContent(iBinX+1,iBinY+1,(histo->GetBinContent(iBinX+1,iBinY+1)+histo_2->GetBinContent(iBinX+1,iBinY+1))/2);
	histo->SetBinError(iBinX+1,iBinY+1,sqrt(histo->GetBinError(iBinX+1,iBinY+1)*histo->GetBinError(iBinX+1,iBinY+1)+
						histo_2->GetBinError(iBinX+1,iBinY+1)*histo_2->GetBinError(iBinX+1,iBinY+1))*0.5);
      }
    }
  }
}


// fix shape uncertainty above xPoint with xValue
void fixShapeUncertainty(TH1* nominalHisto, TH1* sysHisto, float xPoint, float xValue){

  if(sysHisto == 0 || sysHisto == NULL) return;
  for(int iBin = 0; iBin < sysHisto->GetNbinsX(); iBin++){
    if(iBin >= nominalHisto->FindBin(xPoint))
      sysHisto->SetBinContent(iBin+1,nominalHisto->GetBinContent(iBin+1)*xValue);
  }
}

// un-roll a 2D histo into a TH1 with nxm bins from 0 to nxm
// alongX when true means that the 2D histogram is un-rolled looping in x,y order, when false in y,x
TH1* unroll2DHistograms(TH2* nominalHisto, bool alongX = true){

  vector<double> axisLimit;
  for(int binLimit = 0.; binLimit <= nominalHisto->GetNbinsX()*nominalHisto->GetNbinsY(); binLimit++)
    axisLimit.push_back(double(binLimit));

  TH1F* unrolledHist = new TH1F(TString(nominalHisto->GetName()).ReplaceAll("_2D",""),"",axisLimit.size()-1,&axisLimit[0]);
  unrolledHist->Sumw2();
  
  if(alongX){
    int iBin = 1; // for 1D histo
    for(int iBinX = 1; iBinX <= nominalHisto->GetNbinsX(); iBinX++){ // overflow already included by makehist.h
      for(int iBinY = 1; iBinY <= nominalHisto->GetNbinsY(); iBinY++){ 
	unrolledHist->SetBinContent(iBin,nominalHisto->GetBinContent(iBinX,iBinY));
	unrolledHist->SetBinError(iBin,nominalHisto->GetBinError(iBinX,iBinY));
	iBin++;
      }
    }
  }
  else{

    int iBin = 1; // for 1D histo
    for(int iBinY = 1; iBinY <= nominalHisto->GetNbinsY(); iBinY++){ // overflow already included by makehist.h
      for(int iBinX = 1; iBinX <= nominalHisto->GetNbinsX(); iBinX++){ 
	unrolledHist->SetBinContent(iBin,nominalHisto->GetBinContent(iBinX,iBinY));
	unrolledHist->SetBinError(iBin,nominalHisto->GetBinError(iBinX,iBinY));
	iBin++;
      }
    }
  }

  return dynamic_cast<TH1*>(unrolledHist);  
}

// re-construct the 2D histo from the 1D unrolled ones
TH2* roll2DHistograms(TH1* nominalHisto, const string & observable_2D, const int & category, bool alongX = true){

  bin2D bins = selectBinning2D(observable_2D,category);
  TH2F* outputHisto = new TH2F(Form("%s_2D",nominalHisto->GetName()),"",bins.binX.size()-1,&bins.binX[0],bins.binY.size()-1,&bins.binY[0]);

  if(outputHisto->GetNbinsX()*outputHisto->GetNbinsY() != nominalHisto->GetNbinsX())
    cout<<"roll2DHistograms: Huston we have a problem with binning --> please check"<<endl;

  if(alongX){
    int iBinX = 1; 
    for(int iBin = 1; iBin < nominalHisto->GetNbinsX(); iBin++){
      if(iBin <= int((bins.binY.size()-1))){
	outputHisto->SetBinContent(iBinX,iBin,nominalHisto->GetBinContent(iBin));
	outputHisto->SetBinError(iBinX,iBin,nominalHisto->GetBinError(iBin));	
      }
      else{
	outputHisto->SetBinContent(iBinX,iBin-(iBinX-1)*(bins.binY.size()-1),nominalHisto->GetBinContent(iBin));
	outputHisto->SetBinError(iBinX,iBin-(iBinX-1)*(bins.binY.size()-1),nominalHisto->GetBinError(iBin));	
      }
      
      if(iBin % (iBinX*(bins.binY.size()-1)) == 0)
	iBinX++;
    }
  }
  else{

    int iBinY = 1;
    for(int iBin = 1; iBin < nominalHisto->GetNbinsX(); iBin++){
      if(iBin <= int((bins.binX.size()-1))){
        outputHisto->SetBinContent(iBin,iBinY,nominalHisto->GetBinContent(iBin));
        outputHisto->SetBinError(iBin,iBinY,nominalHisto->GetBinError(iBin));
      }
      else{
	outputHisto->SetBinContent(iBin-(iBinY-1)*(bins.binX.size()-1),iBinY,nominalHisto->GetBinContent(iBin));
        outputHisto->SetBinError(iBin-(iBinY-1)*(bins.binX.size()-1),iBinY,nominalHisto->GetBinError(iBin));
      }

      if(iBin > 1 and iBin % (iBinY*(bins.binX.size()-1)) == 0)
        iBinY++;
    }
  }

  return outputHisto;
  
}

// from un-rolled hist, obtain a vector of 1D histograms recovering the real physical binning
// when alongX is false means plot in bin of X, otherwise in bin of Y
vector<TH1F*> transformUnrolledHistogram(TH1* unrolledHisto, const string & observable_2D, const int & category, const bool & alongX = false, const string & suffix = ""){

  vector<TH1F*> outputHistograms;  
  bin2D bins = selectBinning2D(observable_2D,category);

  if(alongX){

    if(unrolledHisto->GetNbinsX() % (bins.binY.size()-1) != 0)
	cerr<<"transforUnrolledHistogram: Problem with binning --> please check"<<endl;
      
    for(size_t iBinX = 1; iBinX <= (bins.binX.size()-1); iBinX++){
      if(suffix == "")
	outputHistograms.push_back(new TH1F(Form("%s_binX_%d",unrolledHisto->GetName(),int(iBinX)),"",bins.binY.size()-1,&bins.binY[0]));
      else
	outputHistograms.push_back(new TH1F(Form("%s_%s_binX_%d",unrolledHisto->GetName(),suffix.c_str(),int(iBinX)),"",bins.binY.size()-1,&bins.binY[0]));
    }

    int nstep = 1;
    for(int iBinX = 1; iBinX <= unrolledHisto->GetNbinsX(); iBinX++){ // loop on x-axis and create a new histogram every bins.binY bins      
      outputHistograms.at(nstep-1)->SetBinContent(iBinX-(nstep-1)*(bins.binY.size()-1),unrolledHisto->GetBinContent(iBinX));
      outputHistograms.at(nstep-1)->SetBinError(iBinX-(nstep-1)*(bins.binY.size()-1),unrolledHisto->GetBinError(iBinX));	

      if(iBinX > 1 and iBinX %(nstep*(bins.binY.size()-1)) == 0)
	nstep++;
      
    }
    return outputHistograms;
  }

  else{ // do the same for bins belonging to the same y-axis one

    if(unrolledHisto->GetNbinsX() % (bins.binX.size()-1) != 0)
      cerr<<"transforUnrolledHistogram: Problem with binning --> please check"<<endl;
      
    for(size_t iBinY = 1; iBinY <= (bins.binY.size()-1); iBinY++){
      if(suffix == "")
	outputHistograms.push_back(new TH1F(Form("%s_binY_%d",unrolledHisto->GetName(),int(iBinY)),"",bins.binX.size()-1,&bins.binX[0]));
      else
	outputHistograms.push_back(new TH1F(Form("%s_%s_binY_%d",unrolledHisto->GetName(),suffix.c_str(),int(iBinY)),"",bins.binX.size()-1,&bins.binX[0]));
    }

    for(size_t iBinY = 1; iBinY <= (bins.binY.size()-1); iBinY++){
      int nstep = iBinY;     
      int ntimes = 1;
      for(int iBinX = 1; iBinX <= unrolledHisto->GetNbinsX(); iBinX++){

	if(iBinX < int(iBinY)) continue;

	if(int(iBinY) == iBinX){
	  outputHistograms.at(iBinY-1)->SetBinContent(nstep-iBinY+1,unrolledHisto->GetBinContent(iBinX));
	  outputHistograms.at(iBinY-1)->SetBinError(nstep-iBinY+1,unrolledHisto->GetBinError(iBinX));
	}
	else if(iBinX == nstep+ntimes*(int(bins.binY.size()-1))){
	  ntimes++;
	  outputHistograms.at(iBinY-1)->SetBinContent(ntimes,unrolledHisto->GetBinContent(iBinX));
	  outputHistograms.at(iBinY-1)->SetBinError(ntimes,unrolledHisto->GetBinError(iBinX));
	}	
      }
    }    
    return outputHistograms;
  }

}

void changeInLatexName(string & variable){

  if(variable == "met")
    variable = "Recoil [GeV]";
  else if(variable == "ht")
    variable = "H_{T} [GeV]";
  else if(variable == "mT")
    variable = "m_{T} [GeV]";
  else if(variable == "njet")
    variable = "N_{jet}";
  else if(variable == "nbjet")
    variable = "N_{bjet}";
  else if(variable == "mpruned")
    variable = "m_{pruned} [GeV]";
  else if(variable == "tau2tau1")
    variable = "#tau_{2}/#tau_{1}";
  else if(variable == "bosonPt")
    variable = "p_{T}^{V} [GeV]";
  else if(variable == "jetPt")
    variable = "p_{T}^{jet} [GeV]";
  else if(variable == "boostedJetPt")
    variable = "p_{T}^{jet} [GeV]";

}


pair<string,string> observableName (string name, bool alongX = false){

  stringstream name_tmp(name.c_str());
  string segment;
  vector<string> seglist;
  while(getline(name_tmp, segment,'_')){
    seglist.push_back(segment);
  }

  string variableX = seglist.back();
  string variableY = seglist.at(seglist.size()-2);

  changeInLatexName(variableX);
  changeInLatexName(variableY);
  

  if(alongX)
    return make_pair(variableX,variableY);
  else
    return make_pair(variableY,variableX);
}


#endif
