#ifndef HISTOUTILS2D_H
#define HISTOUTILS2D_H

#include <vector>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "histoUtils.h"

using namespace std;

////////////////////////////////////////////////////
// define binnings for the different observables // 
///////////////////////////////////////////////////                               

map<string,vector<double> > bins_monoV_2D;
map<string,vector<double> > bins_monoJ_2D;

void initializeBinning2D(){

  bins_monoJ_2D["met"] = {200.,250.,300.,350,400.,500.,600.,750.,950.,1300.};
  bins_monoJ_2D["mpruned"]  = {0.,5.,10.,15.,20.,25.,30.,35.,45.,55.,65.,75.,85.,95.,105.};
  bins_monoJ_2D["tau2tau1"] = {0.,0.15,0.3,0.4,0.5,0.7,0.8,0.9,1.};
  bins_monoJ_2D["njet"] = {1,2,3,10};
  bins_monoJ_2D["njet_v2"] = {1,2,3,4,10};
  bins_monoJ_2D["HT"] = {50.,300.,650.,950.,2000.};
  bins_monoJ_2D["mT"] = {200.,300,400.,500,600.,800,1100,1400,2000.};
  bins_monoJ_2D["dphiJJ"] = {-0.2,0.0,0.6,1.5,3.14};
  bins_monoJ_2D["QGL"] = {0.,0.15,0.50,0.85,1.};
 
}

//// 2D binning
struct bin2D {
  
  vector<double> binX;
  vector<double> binY;

};

bin2D selectBinning2D (const string & observable, const Category & category){

  bin2D bins;
  std::string varX = observable.substr(std::distance(observable.begin(),observable.begin()),observable.find("_"));
  std::string varY = observable.substr(observable.find("_"),std::distance(observable.end(),observable.begin())); 

  if(category == Category::monojet or category == Category::inclusive){
    bins.binX = bins_monoJ_2D[varX];
    bins.binY = bins_monoJ_2D[varY];
    return bins;
  }
  else
    return bins;
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
TH2* roll2DHistograms(TH1* nominalHisto, const string & observable_2D, const Category & category, bool alongX = true){

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
vector<TH1F*> transformUnrolledHistogram(TH1* unrolledHisto, const string & observable_2D, const Category & category, const bool & alongX = false, const string & suffix = ""){

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

/// useful to make interpolated plots TGraph2D correctly
TH2* cloneHistoIncludingOverUnderFlow (TH2* histo){

  vector<double> binX;
  for(int iBinX = 0; iBinX <= histo->GetXaxis()->GetNbins()+2; iBinX++){
    binX.push_back(histo->GetXaxis()->GetBinLowEdge(iBinX));
  }
  vector<double> binY;
  for(int iBinY = 0; iBinY <= histo->GetYaxis()->GetNbins()+2; iBinY++){
    binY.push_back(histo->GetYaxis()->GetBinLowEdge(iBinY));
  }

  TH2F* temp = new TH2F((string(histo->GetName())+"_clone").c_str(),"",int(binX.size()-1),&binX[0],int(binY.size()-1),&binY[0]);
  
  for(int iBinX = 1; iBinX <= temp->GetNbinsX(); iBinX++){
    for(int iBinY = 1; iBinY <= temp->GetNbinsY(); iBinY++){
	temp->SetBinContent(iBinX,iBinY,histo->GetBinContent(histo->FindBin(temp->GetXaxis()->GetBinCenter(iBinX),temp->GetYaxis()->GetBinCenter(iBinY))));
    }
  }
  return dynamic_cast<TH2*>(temp);
}

#endif
