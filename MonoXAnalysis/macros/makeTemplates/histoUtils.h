#ifndef HISTOUTILS_H
#define HISTOUTILS_H

#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TString.h"

using namespace std;

// Basic enum
enum class Sample   { qcd, sig, gam, wmn, zmm, wen, zee, topmu, topel};
enum class Category { inclusive, monojet, monoV, VBF, boosted, prunedMass, tau2tau1};

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

// MET
map<string,vector<double> > bins_monoV; 
map<string,vector<double> > bins_monoJ;
map<string,vector<double> > bins_substructure;
map<string,vector<double> > bins_VBF;

void initializeBinning(){

  bins_monoV["met"] = {250.,300.,350.,400.,500.,600.,750.,1000.};
  bins_monoJ["met"] = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
  bins_substructure["met"] = {200.,225.,250.,275.,300.,325,350.,400.,500.,600.,750.,1000.};
  bins_VBF["met"] = {250.,300.,350.,400.,500.,600.,750.,1000.};

  bins_monoV["jetmetdphi"] = {0,0.25,0.5,0.55,0.65,0.75,0.85,1,1.15,1.30,1.5,1.7,1.9,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.14};
  bins_monoJ["jetmetdphi"] = {0,0.25,0.5,0.55,0.65,0.75,0.85,1,1.15,1.30,1.5,1.7,1.9,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.14};
  bins_VBF["jetmetdphi"]   = {0,0.25,0.5,0.55,0.65,0.75,0.85,1,1.15,1.30,1.5,1.7,1.9,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.14};

  bins_monoV["dphiJJ"] = {-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.14};
  bins_monoJ["dphiJJ"] = {-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.14};
  bins_monoV["minDphiJJ"] = {-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.14};
  bins_monoJ["minDphiJJ"] = {-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.14};
  bins_monoV["minDphiJ1J"] = {-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.14};
  bins_monoJ["minDphiJ1J"] = {-1,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3.14};

  bins_monoV["mT"] = {200.,250.,300.,350.,400.,500.,600.,1000.};
  bins_substructure["mT"] = {200.,250.,300.,350.,400.,500.,600.,1000.};
  bins_monoJ["mT"] = {200.,280,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1350,1550,1750,2000};

  bins_monoV["mpruned"] = {65.,67.5,70.,72.5,75.,77.5,80.,82.5,85.,87.5,90.,92.5,95.,97.5,100.,102.5,105.};
  bins_monoV["mpruned_v2"] = {65.,73.,81.,89.,97.,105.};
  bins_substructure["mpruned"] = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};
  bins_substructure["mpruned_v2"] = {0.,15.,30.,45.,60.,75.,90.,105.,120.,135.,150.,200};
  bins_monoJ["mpruned"] = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};

  bins_monoV["msoftdrop"] = {65.,67.5,70.,72.5,75.,77.5,80.,82.5,85.,87.5,90.,92.5,95.,97.5,100.,102.5,105.};
  bins_monoV["msoftdrop_v2"] = {65.,73.,81.,89.,97.,105.};
  bins_substructure["msoftdrop"] = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};
  bins_substructure["msoftdrop_v2"] = {0.,15.,30.,45.,60.,75.,90.,105.,120.,135.,150.,200};
  bins_monoJ["msoftdrop"] = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};

  bins_monoV["mraw"] = {65.,67.5,70.,72.5,75.,77.5,80.,82.5,85.,87.5,90.,92.5,95.,97.5,100.,102.5,105.};
  bins_monoV["mraw_v2"] = {65.,73.,81.,89.,97.,105.};
  bins_substructure["mraw"] = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};
  bins_substructure["mraw_v2"] = {0.,15.,30.,45.,60.,75.,90.,105.,120.,135.,150.,200};
  bins_monoJ["mraw"] = {0.,7.,14.,21.,28.,35.,42.,49.,56.,63.,70.,77.,84.,91.,98.,105.,112.,119.,126.,133.,140.,147.,154.,161.,168.,175.,182.,189.,196.};

  bins_monoV["njet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_substructure["njet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_monoJ["njet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_VBF["njet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};

  bins_monoV["nbjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_substructure["nbjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_monoJ["nbjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_VBF["nbjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};

  bins_monoV["ncjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_substructure["ncjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_monoJ["ncjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};
  bins_VBF["ncjet"] = {0.,1.,2.,3.,4.,5.,6.,7.,8.};

  bins_monoV["HT"] = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
  bins_monoJ["HT"] = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};
  bins_substructure["HT"] = {150.,200.,250.,300.,350.,400.,450.500,550,600.,650,700.,750,850,950,1050,1250,1450,1650,1850,2100};

  bins_monoV["boostedjetpt"] = {200.,250.,300.,350.,400.,500.,600.,1000.};
  bins_monoV["jetpt"] = {200.,250.,300.,350.,400.,500.,600.,1000.};
  bins_substructure["boostedjetpt"] = {100.,120.,140.,160.,180.,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
  bins_substructure["jetpt"] = {100.,120.,140.,160.,180.,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
  bins_monoJ["jetpt"] = {100.,120.,140.,160.,180.,200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};
  bins_VBF["jetpt"] = {100.,150.,200,250,300,350,400,450,550,650,850};
  bins_VBF["jetpt2"] = {40,60,80,100.,120.,140.,160.,180.,200.,230.,260,350};

  bins_monoV["bosonpt"]    = {200.,250.,300.,350.,400.,500.,600.,1000.}; 
  bins_monoV["bosonpt_v2"] = {50.,70.,90.,120.,150.,180.,210.,230.,250.,300.,350.,400.,500.,600.,1000.}; 
  bins_substructure["bosonpt"] = {50.,70.,90.,120.,150.,180.,210.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250.};
  bins_monoJ["bosonpt"]    = {50.,90.,120.,150.,180.,210.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250.};
  bins_monoJ["bosonpt_v2"] = {200.,250.,300.,350.,400.,500.,600.,1000.};

  bins_monoV["QGL_AK8"] = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};
  bins_substructure["QGL_AK8"] = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};
  bins_monoJ["QGL"] = {0.,0.04,0.08,0.12,0.16,0.24,0.32,0.40,0.48,0.60,0.68,0.76,0.84,0.88,0.92,0.96,1.};

  bins_monoV["tau2tau1"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.90,1.};
  bins_substructure["tau2tau1"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.90,1.};
  bins_monoJ["tau2tau1"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.90,1.};
  
  bins_monoV["nvtx"] = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};
  bins_monoJ["nvtx"] = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};
  bins_substructure["nvtx"] = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};
  bins_VBF["nvtx"] = {0.,2.,4.,6.,8.,10.,12.,14.,16,18,20,22,24,26,28,30,32};

  bins_monoV["chfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_monoJ["chfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_substructure["chfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_VBF["chfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};

  bins_monoV["nhfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_monoJ["nhfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_substructure["nhfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_VBF["nhfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};

  bins_monoV["emfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_monoJ["emfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_substructure["emfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};
  bins_VBF["emfrac"] = {0.,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88,0.92,0.96,1.};

  bins_monoV["btagCSV"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};
  bins_monoJ["btagCSV"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};
  bins_VBF["btagCSV"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};

  bins_monoV["btagCSV_max"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};
  bins_monoJ["btagCSV_max"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};

  bins_monoV["btagCSV_min"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};
  bins_monoJ["btagCSV_min"] = {0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.91,0.94,0.96,0.98,1.};

  bins_VBF["jeteta"]  = {-5,-4.7,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.7,5};
  bins_VBF["jeteta2"] = {-5,-4.7,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.7,5};
  bins_VBF["detajj"]  = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7};
  bins_VBF["mjj"]     = {200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000};

}
 
// binning selections                                                                                                                                                          
vector<double> selectBinning (const string & observable, const Category & category){

  if(category == Category::inclusive or category == Category::monojet)
    return bins_monoJ[observable];
  else if(category == Category::monoV)
    return bins_monoV[observable];
  else if(category == Category::boosted or category == Category::prunedMass or category == Category::tau2tau1)
    return bins_substructure[observable];
  else if(category == Category::VBF)
    return bins_VBF[observable];
  else{
    vector<double> dummy;
    return dummy;
  }
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
TH1F* generateEnvelopeMax(vector<TH1F*> & histoVec, string nameBase){

  if(histoVec.size() == 0){
    cerr<<"generateEnvelopeMax: empty input histogram collection "<<endl;
    TH1F* temp = new TH1F();
    temp->SetName((nameBase+"_envelopeMax").c_str());
    return temp;
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

TH1F* generateEnvelopeMin(vector<TH1F*> & histoVec, string nameBase){

  if(histoVec.size() == 0){
    cerr<<"generateEnvelopeMax: empty input histogram collection "<<endl;
    TH1F* temp = new TH1F();
    temp->SetName((nameBase+"_envelopeMin").c_str());
    return temp;
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
