#include "../makeTemplates/histoUtils.h"

void makeYieldTableFromMLFit(string inputFileName, Category category, bool isZeynep, bool printPrefit = false, bool mergeEWQCD = false){

  TFile* inputFile = TFile::Open(inputFileName.c_str());

  string dir = "shapes_fit_b";
  if(printPrefit) dir = "shapes_prefit";

  TGraphAsymmErrors* datahist = NULL;
  if(isZeynep and category == Category::monojet){
    datahist = (TGraphAsymmErrors*) inputFile->Get((dir+"/monojet_signal/data").c_str());
  }
  else if(not isZeynep and category == Category::monojet){
    datahist = (TGraphAsymmErrors*) inputFile->Get((dir+"/ch1_ch1/data").c_str());
  }
  else if(isZeynep and category == Category::monoV){
    datahist = (TGraphAsymmErrors*) inputFile->Get((dir+"/monov_signal/data").c_str());
  }
  else if(not isZeynep and category == Category::monoV){
    datahist = (TGraphAsymmErrors*) inputFile->Get((dir+"/ch2_ch1/data").c_str());
  }
  else if(category == Category::VBF or category == Category::VBFrelaxed){
    datahist = (TGraphAsymmErrors*) inputFile->Get((dir+"/ch1/data").c_str());
  }

  TH1* zvvhist = NULL;
  TH1* wjethist = NULL;
  TH1* tophist = NULL;
  TH1* dibosonhist = NULL;
  TH1* vghist = NULL;
  TH1* ewkwhist = NULL;
  TH1* ewkzhist = NULL;
  TH1* qcdhist = NULL;
  TH1* dyhist = NULL;
  TH1* gammahist = NULL;
  TH1* totalhist = NULL;

  if(not isZeynep and category == Category::monojet){
    zvvhist  = (TH1*)inputFile->Get((dir+"/ch1_ch1/Znunu").c_str());
    wjethist = (TH1*)inputFile->Get((dir+"/ch1_ch1/WJets").c_str());
    dyhist   = (TH1*)inputFile->Get((dir+"/ch1_ch1/ZJets").c_str());
    tophist  = (TH1*)inputFile->Get((dir+"/ch1_ch1/Top").c_str());
    dibosonhist = (TH1*)inputFile->Get((dir+"/ch1_ch1/Dibosons").c_str());
    vghist   = (TH1*)inputFile->Get((dir+"/ch1_ch1/VGamma").c_str());
    qcdhist = (TH1*)inputFile->Get((dir+"/ch1_ch1/QCD").c_str());
    gammahist = (TH1*)inputFile->Get((dir+"/ch1_ch1/GJets").c_str());
    totalhist = (TH1*)inputFile->Get((dir+"/ch1_ch1/total_background").c_str());
  }
  else if(not isZeynep and category == Category::monoV){
    zvvhist = (TH1*)inputFile->Get((dir+"/ch2_ch1/Znunu").c_str());
    wjethist = (TH1*)inputFile->Get((dir+"/ch2_ch1/WJets").c_str());
    dyhist = (TH1*)inputFile->Get((dir+"/ch2_ch1/ZJets").c_str());
    vghist   = (TH1*)inputFile->Get((dir+"/ch2_ch1/VGamma").c_str());
    tophist = (TH1*)inputFile->Get((dir+"/ch2_ch1/Top").c_str());
    dibosonhist = (TH1*)inputFile->Get((dir+"/ch2_ch1/Dibosons").c_str());
    qcdhist = (TH1*)inputFile->Get((dir+"/ch2_ch1/QCD").c_str());
    gammahist = (TH1*)inputFile->Get((dir+"/ch2_ch1/GJets").c_str());
    totalhist = (TH1*)inputFile->Get((dir+"/ch2_ch1/total_background").c_str());
  }
  else if(isZeynep and category == Category::monojet){
    zvvhist = (TH1*)inputFile->Get((dir+"/monojet_signal/zjets").c_str());
    wjethist = (TH1*)inputFile->Get((dir+"/monojet_signal/wjets").c_str());
    dyhist = (TH1*)inputFile->Get((dir+"/monojet_signal/zll").c_str());
    tophist = (TH1*)inputFile->Get((dir+"/monojet_signal/top").c_str());
    dibosonhist = (TH1*)inputFile->Get((dir+"/monojet_signal/diboson").c_str());
    qcdhist = (TH1*)inputFile->Get((dir+"/monojet_signal/qcd").c_str());
    gammahist = (TH1*)inputFile->Get((dir+"/monojet_signal/gjets").c_str());
    totalhist = (TH1*)inputFile->Get((dir+"/monojet_signal/total_background").c_str());
  }
  else if(isZeynep and category == Category::monoV){
    zvvhist = (TH1*)inputFile->Get((dir+"/monov_signal/zjets").c_str());
    wjethist = (TH1*)inputFile->Get((dir+"/monov_signal/wjets").c_str());
    dyhist = (TH1*)inputFile->Get((dir+"/monov_signal/zll").c_str());
    tophist = (TH1*)inputFile->Get((dir+"/monov_signal/top").c_str());
    dibosonhist = (TH1*)inputFile->Get((dir+"/monov_signal/diboson").c_str());
    qcdhist = (TH1*)inputFile->Get((dir+"/monov_signal/qcd").c_str());
    gammahist = (TH1*)inputFile->Get((dir+"/monov_signal/gjets").c_str());
    totalhist = (TH1*)inputFile->Get((dir+"/monov_signal/total_background").c_str());
  }

  else if(category == Category::VBF or category == Category::VBFrelaxed){
    zvvhist  = (TH1*)inputFile->Get((dir+"/ch1/Znunu").c_str());
    wjethist = (TH1*)inputFile->Get((dir+"/ch1/WJets").c_str());
    dyhist   = (TH1*)inputFile->Get((dir+"/ch1/ZJets").c_str());
    tophist  = (TH1*)inputFile->Get((dir+"/ch1/Top").c_str());
    dibosonhist = (TH1*)inputFile->Get((dir+"/ch1/Dibosons").c_str());
    vghist   = (TH1*)inputFile->Get((dir+"/ch1/VGamma").c_str());
    qcdhist = (TH1*)inputFile->Get((dir+"/ch1/QCD").c_str());
    gammahist = (TH1*)inputFile->Get((dir+"/ch1/GJets").c_str());
    totalhist = (TH1*)inputFile->Get((dir+"/ch1/total_background").c_str());
    ewkzhist = (TH1*)inputFile->Get((dir+"/ch1/Znunu_EWK").c_str()); 
    ewkwhist = (TH1*)inputFile->Get((dir+"/ch1/WJets_EWK").c_str()); 
    
  }

  if(mergeEWQCD){
    if(ewkzhist)
      zvvhist->Add(ewkzhist);
    if(ewkwhist)
      wjethist->Add(ewkwhist);
  }

  TH1* other = (TH1*) qcdhist->Clone("other");
  if(gammahist)
    other->Add(gammahist);
  other->Add(dyhist);

  if(vghist)
    dibosonhist->Add(vghist);

  // make the table
  ofstream outputfile;
  string postfix;
  if(category == Category::monojet)
    postfix = "monojet";
  else if(category == Category::monoV)
    postfix = "monoV";
  else if(category == Category::VBF or category == Category::VBFrelaxed)
    postfix = "VBF";
  
  outputfile.open(Form("yield_%s.txt",postfix.c_str()));
  if(category == Category::monojet or category == Category::monoV or (category == Category::VBF and mergeEWQCD) or (category == Category::VBFrelaxed and mergeEWQCD)){
    outputfile<<"$E_{T}^{miss}$ (GeV) & Observed & $Z \\rightarrow \\nu\\nu$+jets & $W \\rightarrow \\ell\\nu$+jets & Top & Dibosons & Other & Total Bkg. \\\\"<<endl; 
    for(int ibin = 0; ibin < totalhist->GetNbinsX(); ibin++){
      double x,y;
      datahist->GetPoint(ibin,x,y);
      outputfile<<Form("%d-%d",int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)),int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2)))<<" & ";
      outputfile<<Form("%d",int(y*(int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2))-int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)))))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",zvvhist->GetBinContent(ibin+1)*zvvhist->GetBinWidth(ibin+1),zvvhist->GetBinError(ibin+1)*zvvhist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",wjethist->GetBinContent(ibin+1)*wjethist->GetBinWidth(ibin+1),wjethist->GetBinError(ibin+1)*wjethist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",tophist->GetBinContent(ibin+1)*tophist->GetBinWidth(ibin+1),tophist->GetBinError(ibin+1)*tophist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",dibosonhist->GetBinContent(ibin+1)*dibosonhist->GetBinWidth(ibin+1),dibosonhist->GetBinError(ibin+1)*dibosonhist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",other->GetBinContent(ibin+1)*other->GetBinWidth(ibin+1),other->GetBinError(ibin+1)*other->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",totalhist->GetBinContent(ibin+1)*totalhist->GetBinWidth(ibin+1),totalhist->GetBinError(ibin+1)*totalhist->GetBinWidth(ibin+1))<<" \\\\ ";
      outputfile<<"\n";
    }
  }
  else if((category == Category::VBF or category == Category::VBFrelaxed) and not mergeEWQCD){
    outputfile<<"$m_{jj} (GeV) & Observed & Znunu-QCD & Znunu-EWK & W+jets QCD & W+jets EWK & Top & VV+Vgamma & Other & Total Bkg. \\\\ "<<endl;
    for(int ibin = 0; ibin < totalhist->GetNbinsX(); ibin++){
      double x,y;
      datahist->GetPoint(ibin,x,y);
      outputfile<<Form("%d-%d",int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)),int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2)))<<" & ";
      outputfile<<Form("%d",int(y*(int(totalhist->GetXaxis()->GetBinLowEdge(ibin+2))-int(totalhist->GetXaxis()->GetBinLowEdge(ibin+1)))))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",zvvhist->GetBinContent(ibin+1)*zvvhist->GetBinWidth(ibin+1),zvvhist->GetBinError(ibin+1)*zvvhist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",ewkzhist->GetBinContent(ibin+1)*ewkzhist->GetBinWidth(ibin+1),ewkzhist->GetBinError(ibin+1)*ewkzhist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",wjethist->GetBinContent(ibin+1)*wjethist->GetBinWidth(ibin+1),wjethist->GetBinError(ibin+1)*wjethist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",ewkwhist->GetBinContent(ibin+1)*ewkwhist->GetBinWidth(ibin+1),ewkwhist->GetBinError(ibin+1)*ewkwhist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",tophist->GetBinContent(ibin+1)*tophist->GetBinWidth(ibin+1),tophist->GetBinError(ibin+1)*tophist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",dibosonhist->GetBinContent(ibin+1)*dibosonhist->GetBinWidth(ibin+1),dibosonhist->GetBinError(ibin+1)*dibosonhist->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",other->GetBinContent(ibin+1)*other->GetBinWidth(ibin+1),other->GetBinError(ibin+1)*other->GetBinWidth(ibin+1))<<" & ";
      outputfile<<Form("%.3f $\\pm$ %.3f",totalhist->GetBinContent(ibin+1)*totalhist->GetBinWidth(ibin+1),totalhist->GetBinError(ibin+1)*totalhist->GetBinWidth(ibin+1))<<" \\\\ ";
      outputfile<<"\n";      
    }
  }
  outputfile.close();
}

