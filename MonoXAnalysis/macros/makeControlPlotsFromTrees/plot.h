#ifndef PLOT_H
#define PLOT_H

#include "makeGenericHist.h"
#include "drawPlot.h"
#include "datasets.h"

void plot(Sample chan, // channel
	  const char* varstr, // variable
	  const char* cut, // selection
	  const char* xlabel, // label x-axis
	  const char* ylabel, // label y-axis
	  const char* dirname,  // directory name
	  const char* plotname, // plot name
	  double lumi, // lumi
	  int nbins, double* xbins, // binning
	  double ymin, double ymax, 
	  bool isInclusive = false,
	  bool isLog = true) {

    TH1F*  dathist = new TH1F((string(plotname)+"_dat_hist").c_str(), "Data"       , nbins, xbins);
    TH1F*  znnhist = new TH1F((string(plotname)+"_znn_hist").c_str(), "Z(#nu#nu)"  , nbins, xbins);
    TH1F*  wlnhist = new TH1F((string(plotname)+"_wln_hist").c_str(), "W(l#nu)"    , nbins, xbins);
    TH1F*  zllhist = new TH1F((string(plotname)+"_zll_hist").c_str(), "Z(ll)"      , nbins, xbins);
    TH1F*  tophist = new TH1F((string(plotname)+"_top_hist").c_str(), "Top"        , nbins, xbins);
    TH1F*  dibhist = new TH1F((string(plotname)+"_dib_hist").c_str(), "Dibosons"   , nbins, xbins);
    TH1F*  qcdhist = new TH1F((string(plotname)+"_qcd_hist").c_str(), "QCD"        , nbins, xbins);
    TH1F*  gamhist = new TH1F((string(plotname)+"_gam_hist").c_str(), "#gamma+Jets", nbins, xbins);

    double xmin = xbins[0];
    double xmax = xbins[nbins];

    znnhist->SetFillColor(kAzure-9);
    wlnhist->SetFillColor(kYellow-9);
    zllhist->SetFillColor(kOrange-9);
    tophist->SetFillColor(kOrange+2);
    dibhist->SetFillColor(kViolet-9);
    qcdhist->SetFillColor(kGray);
    gamhist->SetFillColor(kCyan-10);

    vector<TH1*> hists;
    hists.push_back(dathist);
    hists.push_back(qcdhist);
    if (chan == Sample::sig || chan == Sample::gam || chan == Sample::topmu || chan == Sample::topel) {
      hists.push_back(gamhist);
    }
    if (chan != Sample::gam and chan != Sample::topmu and chan != Sample::topel) {
      hists.push_back(dibhist);
      hists.push_back(tophist);
    }
    if (chan == Sample::sig || chan == Sample::wmn || chan == Sample::wen || chan == Sample::topmu || chan == Sample::topel) {
      hists.push_back(zllhist);
      hists.push_back(wlnhist);
    }
    if (chan == Sample::zmm || chan == Sample::zee) {
      hists.push_back(wlnhist);
      hists.push_back(zllhist);
    }
    if (chan == Sample::sig) {
      hists.push_back(znnhist);
    }
    if (chan == Sample::topmu || chan == Sample::topel)
      hists.push_back(tophist);
   
    TFile kfactfile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root");
    TH1*  wnlo = (TH1*)kfactfile.Get("EWKcorr/W");
    TH1*  wlo  = (TH1*)kfactfile.Get("WJets_LO/inv_pt");
    TH1*  znlo = (TH1*)kfactfile.Get("EWKcorr/Z");
    TH1*  zlo  = (TH1*)kfactfile.Get("ZJets_LO/inv_pt");
    TH1*  gnlo = (TH1*)kfactfile.Get("EWKcorr/photon");
    TH1*  glo  = (TH1*)kfactfile.Get("GJets_LO/inv_pt_G");
    wnlo->Divide(wlo);
    znlo->Divide(zlo);
    gnlo->Divide(glo);

    vector<TH1*> efac;
    vector<TH1*> zfac;
    vector<TH1*> wfac;
    vector<TH1*> gfac;
    zfac.push_back(znlo);    
    wfac.push_back(wnlo);
    gfac.push_back(gnlo);

    if(isInclusive == true){
      zfac.clear();
      wfac.clear();
      gfac.clear();
    }

    string dtpath     = "/home/rgerosa/MONOJET_ANALYSIS_2016_Data/NoMetCut/";
    string mcpath     = "/home/rgerosa/MONOJET_ANALYSIS_2016_Data/NoMetCut/";
    string dataset    = "";
    string filter     = "";

    if (chan == Sample::sig || chan == Sample::wmn || chan == Sample::zmm || chan == Sample::topmu) dataset = "SingleMuon";
    if (chan == Sample::gam)                                                                        dataset = "SinglePhoton";
    if (chan == Sample::wen || chan == Sample::zee || chan == Sample::topel)                        dataset = "SingleElectron";

    if (chan == Sample::sig) filter = "sig";    
    if (chan == Sample::gam) filter = "gam";    
    if (chan == Sample::wmn) filter = "wmn";    
    if (chan == Sample::zmm) filter = "zmm";    
    if (chan == Sample::wen) filter = "wen";    
    if (chan == Sample::zee) filter = "zee";    
    if (chan == Sample::topmu) filter = "topmu";    
    if (chan == Sample::topel) filter = "topel";    

    vector<string> datfiles;
    vector<string> znnfiles;
    vector<string> wlnfiles;
    vector<string> zllfiles;
    vector<string> topfiles;
    vector<string> dibfiles;
    vector<string> qcdfiles;
    vector<string> gamfiles;

    fillDataSetNames(znnfiles, "znn");
    fillDataSetNames(wlnfiles, "wln");
    fillDataSetNames(zllfiles, "zll");
    fillDataSetNames(topfiles, "top");
    fillDataSetNames(dibfiles, "dib");
    fillDataSetNames(qcdfiles, "qcd");
    fillDataSetNames(gamfiles, "gam");

    datfiles.push_back(dataset);
    datfiles.push_back(dataset);


    for (size_t i = 1; i < datfiles.size(); i++) 
      makeGenericHist(dtpath+"/"+datfiles[0]+"/"+filter+"filter/"+filter+"_"+datfiles[i]+".root" , dathist, varstr, cut, false, chan, lumi, efac, isInclusive);

    if  (chan == Sample::sig) {
      for (size_t i = 1; i < znnfiles.size(); i++) 
	makeGenericHist(mcpath+"/"+znnfiles[0]+"/"+filter+"filter/"+filter+"_"+znnfiles[i]+".root" , znnhist, varstr, cut, true, chan, lumi, zfac, isInclusive);
    }
    if  (chan != Sample::gam) {
      for (size_t i = 1; i < wlnfiles.size(); i++) 
	makeGenericHist(mcpath+"/"+wlnfiles[0]+"/"+filter+"filter/"+filter+"_"+wlnfiles[i]+".root" , wlnhist, varstr, cut,  true, chan, lumi, wfac, isInclusive);
    }
    if  (chan != Sample::gam) {
      for (size_t i = 1; i < zllfiles.size(); i++) 
	makeGenericHist(mcpath+"/"+zllfiles[0]+"/"+filter+"filter/"+filter+"_"+zllfiles[i]+".root" , zllhist, varstr, cut,  true, chan, lumi, zfac, isInclusive);
    }
    if  (chan != Sample::gam) {
      for (size_t i = 1; i < topfiles.size(); i++) 
	makeGenericHist(mcpath+"/"+topfiles[0]+"/"+filter+"filter/"+filter+"_"+topfiles[i]+".root" , tophist, varstr, cut,  true, chan, lumi, efac, isInclusive);
    }
    if  (chan != Sample::gam) {
      for (size_t i = 1; i < dibfiles.size(); i++) 
	makeGenericHist(mcpath+"/"+dibfiles[0]+"/"+filter+"filter/"+filter+"_"+dibfiles[i]+".root" , dibhist, varstr, cut,  true, chan, lumi, efac, isInclusive);
    }
    /*
    if  (chan != Sample::gam) {
    for (size_t i = 1; i < qcdfiles.size(); i++) 
      makeGenericHist(mcpath+"/"+qcdfiles[0]+"/"+filter+"filter/"+filter+"_"+qcdfiles[i]+".root" , qcdhist, varstr, cut,  true, chan, lumi, efac, isInclusive);
    }
    if  (chan == Sample::sig || chan == Sample::gam) {
    for (size_t i = 1; i < gamfiles.size(); i++) 
      makeGenericHist(mcpath+"/"+gamfiles[0]+"/"+filter+"filter/"+filter+"_"+gamfiles[i]+".root" , gamhist, varstr, cut,  true, chan, lumi, gfac, isInclusive);
    }
    */
    if  (chan == Sample::gam) {
    for (size_t i = 1; i < datfiles.size(); i++) 
      makeGenericHist(dtpath+"/"+datfiles[0]+"/"+filter+"filter/"+filter+"_"+datfiles[i]+".root" , qcdhist, varstr, cut, false, Sample::qcd, lumi, efac, isInclusive);
    }
    
    kfactfile.Close();    
    drawPlot(hists, xmin, ymin, xmax, ymax, dirname, plotname, xlabel, ylabel, lumi, true,isLog);
}

void plot(Sample chan, // channel 
	  const char* varstr,  // variable
	  const char* cut,  // selection
	  const char* xlabel, // x-axis
	  const char* ylabel, // y-axis
	  const char* dirname, // directory for plot
	  const char* plotname, // plot name 
	  double lumi,  // lumi
	  int    nbins, double xmin, double xmax, // binning
	  double ymin, double ymax, // y-axis
	  bool   isInclusive = false,// run inclusive analysis
	  bool   isLog = true
	  ) {

    double* xbins = new double[nbins+1];
    double binwidth = (xmax - xmin) / double(nbins);
    for (int i = 0; i <= nbins; i++) xbins[i] = xmin + double(i)*binwidth;        
    plot(chan, varstr, cut, xlabel, ylabel, dirname, plotname, lumi, nbins, xbins, ymin, ymax,isInclusive,isLog);

}

#endif
