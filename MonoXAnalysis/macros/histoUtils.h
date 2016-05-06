#ifndef HISTOUTILS_H
#define HISTOUTILS_H

#include <vector>
#include <fstream>
#include <sstream>
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

////////////////////////////////////////////////////
// define binnings for the different observables // 
///////////////////////////////////////////////////                               

vector<double> bins_monoV_met         = {250.,300.,350.,400.,500.,600.,750.,1000.};
//vector<double> bins_monoV_met         = {250.,260.,275.,300.,325.,350.,375.,400.,450.,500.,550.,600.,700.,800.,1000.};
vector<double> bins_monoV_met_v2      = {250.,300.,350.,400.,500.,600.,700.,1000.};
vector<double> bins_monoV_met_v3      = {250.,300.,350.,400.,500.,600.,800.,1000.};
vector<double> bins_substructure_met  = {200.,225.,250.,275.,300.,325,350.,400.,500.,600.,750.,1000.};
vector<double> bins_monoJ_met         = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
vector<double> bins_monoJ_met_v2      = {200.,250.,300.,350.,400.,500.,600.,1000.};

vector<double> bins_monoJ_dphiJJ      = {-0.1,0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.14};
vector<double> bins_monoV_dphiJJ      = {-0.1,0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.14};

vector<double> bins_monoV_mT          = {200.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_substructure_mT   = {200.,250.,300.,350.,400.,500.,600.,1000.};
vector<double> bins_monoJ_mT          = {200.,280,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1350,1550,1750,2000};

vector<double> bins_monoV_mpr         = {65.,67.5,70.,72.5,75.,77.5,80.,82.5,85.,87.5,90.,92.5,95.,97.5,100.,102.5,105.};
vector<double> bins_monoV_mpr_v2      = {65.,73.,81.,89.,97.,105.};

vector<double> bins_monoJ_mpr         = {0.,3.,6.,9.,12.,15.,18.,21.,24.,27.,30.,33.,36.,39.,42.,45.,48.,51.,54.,57.,60.,64.,68.,72.,76.,80.,84.,88.,92.,96.,100.};
vector<double> bins_substructure_mpr  = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};
vector<double> bins_substructure_mpr_v2  = {0.,15.,30.,45.,60.,75.,90.,105.,120.,135.,150.,200};

vector<double> bins_monoV_njet        = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
vector<double> bins_monoJ_njet        = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
vector<double> bins_substructure_njet = {0.,1.,2.,3.,4.,5.,6.,7.,8.};

vector<double> bins_monoV_HT          = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
vector<double> bins_monoJ_HT          = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
vector<double> bins_substructure_HT   = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};

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
  else if((observable == "mpruned_v2" or observable == "msoftdrop_v2" or observable == "mraw") and category >= 3)
    return bins_substructure_mpr_v2;
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
  
  else if((observable == "QGL" or TString(observable).Contains("QGL_")) and category <= 1)
    return bins_monoJ_QGL;
  else if((observable == "QGL" or TString(observable).Contains("QGL_")) and category > 1 and category <= 3)
    return bins_monoV_QGL;
  else if((observable == "QGL" or TString(observable).Contains("QGL_")) and category > 3)
    return bins_monoV_QGL;

  else if((observable == "dphiJJ" or  observable == "minDphiJJ" or observable == "minDphiJ1J") and category <= 1)
    return bins_monoJ_dphiJJ;
  else if((observable == "dphiJJ" or  observable == "minDphiJJ" or observable == "minDphiJ1J") and category > 1)
    return bins_monoV_dphiJJ;
  

  vector<double> dummy;
  return dummy;

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

// make the average of two histograms (TH1)
void makeAverage(TH1* histo, TH1* histo_2){

  if(histo_2->Integral() != 0){ // sanity check to understand if the second has been filled
    for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
      histo->SetBinContent(iBin+1,(histo->GetBinContent(iBin+1)+histo_2->GetBinContent(iBin+1))/2);
      histo->SetBinError(iBin+1,sqrt(histo->GetBinError(iBin+1)*histo->GetBinError(iBin+1)+histo_2->GetBinError(iBin+1)*histo_2->GetBinError(iBin+1))*0.5);
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

void fixShapeUncertainty(TH1* nominalHisto, TH1* sysHisto, int xBin, float xValue){

  if(sysHisto == 0 || sysHisto == NULL) return;
  for(int iBin = 0; iBin < sysHisto->GetNbinsX(); iBin++){
    if(iBin >= xBin)
      sysHisto->SetBinContent(iBin+1,nominalHisto->GetBinContent(iBin+1)*xValue);
  }
}

// dummy bin content
void addDummyBinContent(TH1* histo){

  cerr<<"addDummyBinContent: called for histo "<<histo->GetName()<<endl;
  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    histo->SetBinContent(iBin+1,10.e-6);
    histo->SetBinError(iBin+1,10.e-6);
  }

}

// maxium envelop between histograms
TH1F* generateEnvelopeMax(vector<TH1F*> & histoVec){

  if(histoVec.size() == 0){
    cerr<<"generateEnvelopeMax: empty input histogram collection "<<endl;
    return new TH1F();
  }

  TH1F* histo = (TH1F*) histoVec.at(0)->Clone("envelopeMax");
  histo->Reset();

  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    for(auto hist : histoVec){
      if(hist->GetBinContent(iBin+1) > histo->GetBinContent(iBin+1) and histo->GetBinContent(iBin+1) != 0){
        histo->SetBinContent(iBin+1,hist->GetBinContent(iBin+1));
        histo->SetBinError(iBin+1,hist->GetBinError(iBin+1));
      }
      else if(histo->GetBinContent(iBin+1) == 0){
        histo->SetBinContent(iBin+1,hist->GetBinContent(iBin+1));
        histo->SetBinError(iBin+1,hist->GetBinError(iBin+1));
      }
    }
  }

  return histo;

}

TH1F* generateEnvelopeMin(vector<TH1F*> & histoVec){

  if(histoVec.size() == 0){
    cerr<<"generateEnvelopeMax: empty input histogram collection "<<endl;
    return new TH1F();
  }

  TH1F* histo = (TH1F*) histoVec.at(0)->Clone("envelopeMin");
  histo->Reset();

  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    for(auto hist : histoVec){
      if(hist->GetBinContent(iBin+1) < histo->GetBinContent(iBin+1) and histo->GetBinContent(iBin+1) != 0){
        histo->SetBinContent(iBin+1,hist->GetBinContent(iBin+1));
        histo->SetBinError(iBin+1,hist->GetBinError(iBin+1));
      }
      else if(histo->GetBinContent(iBin+1) == 0){
        histo->SetBinContent(iBin+1,hist->GetBinContent(iBin+1));
        histo->SetBinError(iBin+1,hist->GetBinError(iBin+1));
      }
    }
  }

  return histo;

}

#endif
